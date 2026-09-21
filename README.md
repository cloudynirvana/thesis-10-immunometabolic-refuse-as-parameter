# Immunometabolic tumour-immune interaction ODEs under explicit non-parameters: lactate, checkpoint proxies, and host constraints that must not enter Θ

**Thesis #10 (was R3).** Computational research thesis, in Nile University B.Sc. chapter form, for journal or thesis handoff.

**Author:** Kelechi Emeka Ogbonna  
**Email:** kelechiogbonna300@gmail.com  
**Affiliation:** Independent computational research  
**Date:** 21 September 2026

The chapter headings follow a B.Sc. project skeleton. The calculations are an in-silico ledger on a toy ODE. They are not a wet-lab study and not a clinical series.

## Non-claims

This thesis is research only. It is not a medical device, not clinical decision support, not a protocol, not a dose, and not a cure. It is not checkpoint dosing advice. Lactate, the checkpoint proxy, and the host bound are evidence objects. They are not kinetic parameters and not prescriptions.

No DOI is registered for this document. Do not invent one.

See [DISCLAIMER.md](DISCLAIMER.md).

## The problem (one sentence)

Can lactate, checkpoint proxies, and host constraints stay as refuse-as-parameter evidence objects that change hypothesis rank without entering Θ as kinetic parameters?

## Files

| Path | Role |
| --- | --- |
| `THESIS.md` | Full manuscript (Vancouver citations) |
| `THESIS.pdf` | PDF built from the Markdown |
| `build_pdf.py` | Regenerates `THESIS.pdf` |
| `CITATION.cff` | Citation metadata (no document DOI) |
| `DISCLAIMER.md` | Research-only disclaimer |
| `sim/refuse_gate.py` | Toy ODE, hypothesis rank, and refusal gate (seeds 20260921 to 20260923) |
| `sim/results.json` | Numbers reported in Chapter Four |
| `sim/ledger.json` | Evidence objects and refused proposals |
| `sim/figures/` | Trajectory, rank, and residual figures |

## Reproduce the sketch

```bash
python3 -m pip install -r sim/requirements.txt
python3 sim/refuse_gate.py
python3 build_pdf.py
```

NumPy, SciPy, and Matplotlib are required for the sketch. WeasyPrint and the Python `markdown` package are required only to rebuild the PDF. Trajectories are computational signatures, not patient outcomes.

## How to cite

Ogbonna KE. Immunometabolic tumour-immune interaction ODEs under explicit non-parameters: lactate, checkpoint proxies, and host constraints that must not enter Θ [Internet]. Thesis #10 computational research thesis. 21 September 2026 [cited YYYY Mon DD]. Available from: https://github.com/cloudynirvana/thesis-10-immunometabolic-refuse-as-parameter

Prefer `CITATION.cff` for machine-readable citation. If a document DOI is minted later, add it there only after it exists.

Catalogue: [research-theses-hub](https://github.com/cloudynirvana/research-theses-hub).

## Licence

Manuscript text and sketch code are released under MIT. Computational research only. Not a care product.
