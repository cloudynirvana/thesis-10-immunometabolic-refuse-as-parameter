#!/usr/bin/env python3
"""Toy tumour-immune-lactate ODE and a refuse-as-parameter gate.

Research sketch only. States and coefficients are dimensionless cartoons.
They are not concentrations, doses, schedules, or patient outcomes.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parent
FIG = ROOT / "figures"
SEED_BASE = 20260921

# Kinetic vector Theta. Ledger constants are intentionally absent.
THETA_NAMES = ("r", "K", "kappa", "sigma", "delta", "pi", "lam")
REFUSED_NAMES = (
    "beta",
    "gamma",
    "eta",
    "phi",
    "L_lo",
    "L_hi",
    "sigma_host",
    "checkpoint",
)

THETA_TRUE = {
    "r": 0.30,
    "K": 1.20,
    "kappa": 0.50,
    "sigma": 0.12,
    "delta": 0.25,
    "pi": 0.70,
    "lam": 0.55,
}
X0 = np.array([0.25, 0.30, 0.05], dtype=float)

# Declared on the ledger before any ranking. Not elements of Theta.
PHI = 0.55
SIGMA_HOST = 0.20
SIGMA_SUPPLY_ALT = 0.35
BAND_LOOSE = (0.00, 1.50)
BAND_TIGHT = (0.00, 0.45)

# Synthetic generator used only to tempt promotion. beta is not accepted into Theta.
BETA_GENERATOR = 0.35

T_END = 40.0
N_OBS = 21
CV = 0.03
OBS_FLOOR = 0.05
TIMES = np.linspace(0.0, T_END, N_OBS)


def rhs(t, x, theta, structure, ledger):
    T, E, L = x
    r = theta["r"]
    K = theta["K"]
    kappa = theta["kappa"]
    sigma = theta["sigma"]
    delta = theta["delta"]
    pi = theta["pi"]
    lam = theta["lam"]
    kill = kappa
    if structure == "M_ck":
        kill = ledger["phi"] * kappa
    growth = r * T * (1.0 - T / K) - kill * E * T
    supply = sigma
    if structure == "M_illegal":
        checkpoint = ledger.get("checkpoint_numeric", 0.0)
        gamma = theta.get("gamma", 0.0)
        beta = theta.get("beta", 0.0)
        supply = sigma * (1.0 - gamma * checkpoint)
        effector = supply - delta * E - beta * L * E
    else:
        effector = supply - delta * E
    lactate = pi * T - lam * L
    return (growth, effector, lactate)


def simulate(theta, structure, ledger, t_eval=TIMES):
    sol = solve_ivp(
        lambda t, x: rhs(t, x, theta, structure, ledger),
        (float(t_eval[0]), float(t_eval[-1])),
        X0,
        t_eval=t_eval,
        rtol=1e-7,
        atol=1e-9,
        dense_output=False,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y.T


def observe(traj, rng):
    noise = rng.normal(0.0, 1.0, size=traj.shape)
    return traj * (1.0 + CV * noise)


def rss(y_obs, y_hat):
    scale = np.maximum(np.abs(y_obs), OBS_FLOOR)
    resid = (y_obs - y_hat) / scale
    return float(np.sum(resid ** 2))


def terminal_lactate(theta, structure, ledger):
    return float(simulate(theta, structure, ledger)[-1, 2])


def admissible(theta, structure, ledger, band):
    reasons = []
    if theta["sigma"] > ledger["sigma_host"] + 1e-12:
        reasons.append("sigma_above_host_budget")
    L_end = terminal_lactate(theta, structure, ledger)
    lo, hi = band
    if L_end < lo - 1e-9 or L_end > hi + 1e-9:
        reasons.append("terminal_lactate_outside_host_band")
    return len(reasons) == 0, reasons, L_end


def propose_theta(names, ledger_doc):
    """Return False when a proposal would write a refused symbol into Theta."""
    hit = [n for n in names if n in REFUSED_NAMES]
    record = {
        "proposed_names": list(names),
        "refused_hits": hit,
        "accepted": len(hit) == 0,
        "written_theta": list(THETA_NAMES) if hit else list(names),
    }
    ledger_doc["proposals"].append(record)
    return record


def rank_scenario(y_obs, candidates, ledger, band):
    rows = []
    for name, theta, structure in candidates:
        ok, reasons, L_end = admissible(theta, structure, ledger, band)
        y_hat = simulate(theta, structure, ledger)
        score = rss(y_obs, y_hat)
        rows.append(
            {
                "hypothesis": name,
                "structure": structure,
                "admissible": ok,
                "reasons": reasons,
                "terminal_L": L_end,
                "rss": score,
                "theta_names": list(theta.keys()),
                "theta_length": len([k for k in theta if k in THETA_NAMES or k in ("beta", "gamma", "eta")]),
            }
        )
    admissible_rows = [r for r in rows if r["admissible"]]
    admissible_rows.sort(key=lambda r: r["rss"])
    for i, r in enumerate(admissible_rows, start=1):
        r["rank"] = i
    for r in rows:
        if not r["admissible"]:
            r["rank"] = None
    return rows


def pack_theta(values, names):
    return {n: float(v) for n, v in zip(names, values)}


def fit_illegal(y_obs, ledger):
    """Least squares on Theta plus beta. gamma is held out because C = 0 in this generator."""
    names = list(THETA_NAMES) + ["beta"]
    x0 = np.array([THETA_TRUE[n] for n in THETA_NAMES] + [0.0], dtype=float)
    lower = np.array([0.05, 0.4, 0.05, 0.02, 0.05, 0.1, 0.1, 0.0])
    upper = np.array([1.0, 3.0, 2.0, 0.5, 1.0, 2.0, 2.0, 2.0])

    def fun(v):
        theta = pack_theta(v, names)
        try:
            y_hat = simulate(theta, "M_illegal", ledger)
        except Exception:
            return np.ones(y_obs.size) * 1e3
        scale = np.maximum(np.abs(y_obs), OBS_FLOOR)
        return ((y_obs - y_hat) / scale).ravel()

    fit = least_squares(fun, x0, bounds=(lower, upper), max_nfev=40, ftol=1e-8)
    theta_hat = pack_theta(fit.x, names)
    y_hat = simulate(theta_hat, "M_illegal", ledger)
    return {
        "success": bool(fit.success),
        "nfev": int(fit.nfev),
        "cost": float(fit.cost),
        "rss": rss(y_obs, y_hat),
        "theta_hat": theta_hat,
        "message": str(fit.message),
    }


def fisher_legal(ledger):
    """Numerical Fisher matrix for the legal parameter vector at THETA_TRUE."""
    y = simulate(THETA_TRUE, "M_base", ledger)
    eps = 1e-4
    cols = []
    base_vec = np.array([THETA_TRUE[n] for n in THETA_NAMES])
    for i, name in enumerate(THETA_NAMES):
        step = eps * max(abs(base_vec[i]), 1e-2)
        up = dict(THETA_TRUE)
        dn = dict(THETA_TRUE)
        up[name] = base_vec[i] + step
        dn[name] = base_vec[i] - step
        yup = simulate(up, "M_base", ledger)
        ydn = simulate(dn, "M_base", ledger)
        cols.append(((yup - ydn) / (2.0 * step)).ravel())
    J = np.column_stack(cols)
    scale = np.maximum(np.abs(y.ravel()), OBS_FLOOR)
    var = (CV * scale) ** 2
    W = 1.0 / var
    F = J.T @ (W[:, None] * J)
    evals, evecs = np.linalg.eigh(F)
    evals = np.clip(evals, 0.0, None)
    order = np.argsort(evals)[::-1]
    evals_desc = evals[order]
    weak = evecs[:, order[-1]]
    weak = weak / np.linalg.norm(weak)
    tol = 1e-6 * evals_desc[0]
    rank = int(np.sum(evals_desc > tol))

    # Illegal column: sensitivity to beta at beta = 0, then the gate drops it.
    beta_theta = dict(THETA_TRUE)
    beta_theta["beta"] = 0.0
    beta_theta["gamma"] = 0.0
    ledger_b = dict(ledger)
    ledger_b["checkpoint_numeric"] = 0.0
    up = dict(beta_theta)
    up["beta"] = 1e-3
    y0 = simulate(beta_theta, "M_illegal", ledger_b)
    y1 = simulate(up, "M_illegal", ledger_b)
    j_beta = ((y1 - y0) / 1e-3).ravel()
    J_ext = np.column_stack([J, j_beta])
    F_ext = J_ext.T @ (W[:, None] * J_ext)
    ev_ext = np.sort(np.clip(np.linalg.eigvalsh(F_ext), 0.0, None))[::-1]
    rank_ext = int(np.sum(ev_ext > 1e-6 * ev_ext[0]))
    # Drop the refused column and recompute. This is the gate, not a new model.
    J_kept = J_ext[:, : len(THETA_NAMES)]
    F_kept = J_kept.T @ (W[:, None] * J_kept)
    ev_kept = np.sort(np.clip(np.linalg.eigvalsh(F_kept), 0.0, None))[::-1]
    rank_kept = int(np.sum(ev_kept > 1e-6 * ev_kept[0]))
    return {
        "parameter_order": list(THETA_NAMES),
        "eigenvalues_desc": [float(v) for v in evals_desc],
        "weak_eigenvector": {n: float(c) for n, c in zip(THETA_NAMES, weak)},
        "rank": rank,
        "dimension": len(THETA_NAMES),
        "condition_number": float(evals_desc[0] / max(evals_desc[-1], 1e-30)),
        "illegal_dimension_if_beta_added": int(J_ext.shape[1]),
        "rank_if_beta_column_kept": rank_ext,
        "rank_after_gate_drops_beta_column": rank_kept,
        "tolerance_rule": "eigenvalue > 1e-6 * largest eigenvalue",
        "tolerance_value": float(tol),
    }


def style_ax(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=8)
    ax.grid(True, axis="y", color="#dddddd", linewidth=0.6)


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    ledger_constants = {
        "phi": PHI,
        "sigma_host": SIGMA_HOST,
        "sigma_supply_alt": SIGMA_SUPPLY_ALT,
        "band_loose": list(BAND_LOOSE),
        "band_tight": list(BAND_TIGHT),
        "checkpoint_numeric": 0.0,
    }
    ledger_doc = {
        "object": "refuse-as-parameter ledger",
        "refused_names": list(REFUSED_NAMES),
        "legal_theta_names": list(THETA_NAMES),
        "constants": {
            "phi": PHI,
            "sigma_host": SIGMA_HOST,
            "sigma_supply_alt": SIGMA_SUPPLY_ALT,
            "band_loose": list(BAND_LOOSE),
            "band_tight": list(BAND_TIGHT),
        },
        "evidence_objects": [],
        "proposals": [],
    }

    y_base = simulate(THETA_TRUE, "M_base", ledger_constants)
    y_ck = simulate(THETA_TRUE, "M_ck", ledger_constants)

    theta_supply = dict(THETA_TRUE)
    theta_supply["sigma"] = SIGMA_SUPPLY_ALT

    theta_gen = dict(THETA_TRUE)
    theta_gen["beta"] = BETA_GENERATOR
    theta_gen["gamma"] = 0.0
    y_gen = simulate(theta_gen, "M_illegal", ledger_constants)

    rng0 = np.random.default_rng(SEED_BASE)
    rng1 = np.random.default_rng(SEED_BASE + 1)
    rng2 = np.random.default_rng(SEED_BASE + 2)
    d0 = observe(y_base, rng0)
    d1 = observe(y_ck, rng1)
    d2 = observe(y_gen, rng2)

    ledger_doc["evidence_objects"] = [
        {
            "id": "ev.lactate",
            "kind": "trajectory",
            "channels": ["L"],
            "role": "recorded series and host band; not a kinetic coefficient",
            "bands": {"loose": list(BAND_LOOSE), "tight": list(BAND_TIGHT)},
        },
        {
            "id": "ev.checkpoint",
            "kind": "categorical_proxy",
            "levels": ["absent", "present"],
            "role": "selects whether M_ck is in the candidate set; not a rate and not a dose",
            "phi_if_present": PHI,
        },
        {
            "id": "ev.host",
            "kind": "inequality",
            "statement": "sigma <= sigma_host",
            "sigma_host": SIGMA_HOST,
            "role": "feasibility mask; the ceiling is not estimated",
        },
    ]

    candidates_absent = [
        ("M_base", THETA_TRUE, "M_base"),
        ("M_supply", theta_supply, "M_base"),
    ]
    candidates_present = [
        ("M_base", THETA_TRUE, "M_base"),
        ("M_ck", THETA_TRUE, "M_ck"),
        ("M_supply", theta_supply, "M_base"),
    ]

    scenarios = {
        "checkpoint_absent_band_loose": rank_scenario(d0, candidates_absent, ledger_constants, BAND_LOOSE),
        "checkpoint_present_band_loose": rank_scenario(d1, candidates_present, ledger_constants, BAND_LOOSE),
        "checkpoint_present_band_tight": rank_scenario(d1, candidates_present, ledger_constants, BAND_TIGHT),
    }

    # Illegal promotion attempts. None of these writes are accepted.
    blocked_beta = propose_theta(list(THETA_NAMES) + ["beta", "gamma", "eta"], ledger_doc)
    blocked_phi = propose_theta(list(THETA_NAMES) + ["phi"], ledger_doc)
    blocked_host = propose_theta(list(THETA_NAMES) + ["sigma_host", "L_hi"], ledger_doc)
    allowed = propose_theta(list(THETA_NAMES), ledger_doc)

    legal_on_temptation = rss(d2, y_base)
    generator_on_temptation = rss(d2, y_gen)
    # Frozen M_ck scored on the temptation data, for comparison only.
    mck_on_temptation = rss(d2, y_ck)
    fit = fit_illegal(d2, ledger_constants)
    # The fit is a calculation. It is not written into the accepted vector.
    fit_proposal = propose_theta(list(THETA_NAMES) + ["beta"], ledger_doc)

    fisher = fisher_legal(ledger_constants)

    functionals = {
        "t_grid": [float(t) for t in TIMES],
        "M_base_terminal": {
            "T": float(y_base[-1, 0]),
            "E": float(y_base[-1, 1]),
            "L": float(y_base[-1, 2]),
        },
        "M_ck_terminal": {
            "T": float(y_ck[-1, 0]),
            "E": float(y_ck[-1, 1]),
            "L": float(y_ck[-1, 2]),
        },
        "generator_terminal": {
            "T": float(y_gen[-1, 0]),
            "E": float(y_gen[-1, 1]),
            "L": float(y_gen[-1, 2]),
            "beta": BETA_GENERATOR,
        },
        "rss_temptation": {
            "legal_M_base_frozen": legal_on_temptation,
            "legal_M_ck_frozen": mck_on_temptation,
            "generator_with_beta": generator_on_temptation,
            "least_squares_theta_plus_beta": fit["rss"],
        },
    }

    results = {
        "seed_base": SEED_BASE,
        "seeds": {"D_absent": SEED_BASE, "D_present": SEED_BASE + 1, "D_temptation": SEED_BASE + 2},
        "cv": CV,
        "obs_floor": OBS_FLOOR,
        "t_end": T_END,
        "n_obs": N_OBS,
        "x0": [float(v) for v in X0],
        "theta_true": THETA_TRUE,
        "theta_length_legal": len(THETA_NAMES),
        "functionals": functionals,
        "scenarios": scenarios,
        "gate": {
            "blocked_beta_gamma_eta": blocked_beta,
            "blocked_phi": blocked_phi,
            "blocked_host_symbols": blocked_host,
            "allowed_legal": allowed,
            "blocked_fit_write": fit_proposal,
        },
        "illegal_fit": fit,
        "fisher": fisher,
        "ledger": ledger_doc,
        "units": "dimensionless cartoon; time in arbitrary days",
    }

    (ROOT / "results.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    (ROOT / "ledger.json").write_text(json.dumps(ledger_doc, indent=2), encoding="utf-8")

    # Figure 1. Trajectories.
    fig, axes = plt.subplots(1, 3, figsize=(9.2, 3.15), sharex=True)
    labels = [("T", 0, "Tumour burden"), ("E", 1, "Effector state"), ("L", 2, "Lactate state")]
    for ax, (key, idx, title) in zip(axes, labels):
        ax.plot(TIMES, y_base[:, idx], color="#1b4f72", lw=1.8, label="M_base")
        ax.plot(TIMES, y_ck[:, idx], color="#b03a2e", lw=1.8, ls="--", label="M_ck")
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("Time (cartoon days)", fontsize=8)
        style_ax(ax)
    axes[0].set_ylabel("State (dimensionless)", fontsize=8)
    axes[0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "trajectories_states.png", dpi=150)
    plt.close(fig)

    # Figure 2. Rank scores. Inadmissible bars are hatched.
    fig, axes = plt.subplots(1, 3, figsize=(9.2, 3.4), sharey=True)
    scenario_order = [
        ("checkpoint_absent_band_loose", "Checkpoint absent\nloose lactate band"),
        ("checkpoint_present_band_loose", "Checkpoint present\nloose lactate band"),
        ("checkpoint_present_band_tight", "Checkpoint present\ntight lactate band"),
    ]
    for ax, (key, title) in zip(axes, scenario_order):
        rows = scenarios[key]
        names = [r["hypothesis"] for r in rows]
        vals = [r["rss"] for r in rows]
        colors = ["#1b4f72" if r["admissible"] else "#b0b0b0" for r in rows]
        bars = ax.bar(names, vals, color=colors, width=0.72)
        for bar, r in zip(bars, rows):
            if not r["admissible"]:
                bar.set_hatch("//")
            tag = f"rank {r['rank']}" if r["rank"] else "out"
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), tag, ha="center", va="bottom", fontsize=7)
        ax.set_title(title, fontsize=9)
        style_ax(ax)
    axes[0].set_ylabel("Scaled residual sum of squares", fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "hypothesis_rank.png", dpi=150)
    plt.close(fig)

    # Figure 3. Temptation versus refusal.
    fig, ax = plt.subplots(figsize=(6.6, 3.5))
    names = ["M_base\nadmitted", "M_ck\nnot in set", "Generator\nrefused", "Fit + beta\nrefused"]
    vals = [
        legal_on_temptation,
        mck_on_temptation,
        generator_on_temptation,
        fit["rss"],
    ]
    colors = ["#1b4f72", "#b0b0b0", "#922b21", "#922b21"]
    bars = ax.bar(names, vals, color=colors, width=0.72)
    bars[1].set_hatch("//")
    bars[2].set_hatch("xx")
    bars[3].set_hatch("xx")
    ax.set_ylabel("Scaled RSS on the temptation series", fontsize=8)
    ax.set_title("A lower residual is not membership in Theta", fontsize=10)
    style_ax(ax)
    fig.tight_layout()
    fig.savefig(FIG / "promotion_rss.png", dpi=150)
    plt.close(fig)

    winner = {}
    for key, rows in scenarios.items():
        ranked = [r for r in rows if r["rank"] == 1]
        winner[key] = ranked[0]["hypothesis"] if ranked else None
    summary = {
        "winners": winner,
        "terminals": functionals,
        "fit_beta": fit["theta_hat"].get("beta"),
        "fit_rss": fit["rss"],
        "legal_rss": legal_on_temptation,
        "generator_rss": generator_on_temptation,
        "fisher_rank": fisher["rank"],
        "gate_accepted_flags": [p["accepted"] for p in ledger_doc["proposals"]],
    }
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
