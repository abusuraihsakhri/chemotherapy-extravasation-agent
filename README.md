# Chemotherapy Extravasation Reference Tool

### [Open the Live Application →](https://abusuraihsakhri.github.io/chemotherapy-extravasation-agent/)

A Python and browser reference implementation for antineoplastic extravasation assessment. It provides drug classification, CTCAE v5.0 infusion-site extravasation mapping, agent-specific thermal guidance, selected antidote calculations, an ordered response checklist, and CSV batch processing.

> **Clinical use:** Verify the antineoplastic agent, concentration, institutional extravasation protocol, current product label, and pharmacy preparation before treatment. The software is a reference aid and is not a substitute for urgent clinical assessment or local policy.

## Features

- Curated registry of 25 antineoplastic agents with vesicant/irritant classification and thermal-management metadata.
- CTCAE v5.0 **Infusion site extravasation** mapping using the published term criteria rather than numeric pain cutoffs.
- Dexrazoxane dosing with BSA calculation, dose caps, six-hour timing check, and 50% reduction when creatinine clearance is below 40 mL/min.
- Hyaluronidase, sodium thiosulfate, and topical DMSO reference calculations with agent-specific guardrails.
- Explicit cisplatin concentration/volume threshold handling; missing concentration is not inferred.
- No assumed adult BSA or DMSO treatment area when measurements are missing.
- CLI assessment, isolated antidote calculation, CTCAE mapping, interactive use, and CSV batch processing.
- Browser UI using the same Python engine through Pyodide/WebAssembly, with light/dark themes and no server-side clinical processing.

The legacy `risk` CLI command is retained for compatibility. Its 0–100 score is an **unvalidated heuristic**, not a clinical prediction rule, and should not determine treatment or vascular-access decisions.

## Browser use

The static application is defined by `index.html`, `styles.css`, and `app.js`. On first load it downloads the pinned Pyodide runtime, then executes `chemotherapy_extravasation.py` locally in the browser.

The browser form deliberately does not request a patient identifier. Event inputs and results remain in browser memory and are not sent to an application backend. Only the color-theme preference is saved in local storage. Loading Pyodide requires network access to jsDelivr.

Modern browsers with WebAssembly support are required.

## Command line

Python 3.9 or later is supported.

```bash
python -m pip install .
chemotherapy-extravasation --help
```

Example assessment:

```bash
chemotherapy-extravasation assess \
  --drug doxorubicin \
  --catheter peripheral_forearm \
  --volume 8 \
  --elapsed 1.2 \
  --height 172 \
  --weight 68.5 \
  --age 58 \
  --creatinine 1.0 \
  --pain 5 \
  --json
```

For cisplatin, provide concentration when evaluating the high-concentration pathway:

```bash
chemotherapy-extravasation assess \
  --drug cisplatin \
  --catheter peripheral_forearm \
  --volume 25 \
  --concentration 0.6 \
  --elapsed 0.5
```

Dexrazoxane calculation requires either a verified BSA or both height and weight:

```bash
chemotherapy-extravasation antidote \
  --type dexrazoxane \
  --height 175 \
  --weight 70 \
  --crcl 55 \
  --elapsed 1
```

Batch processing:

```bash
chemotherapy-extravasation batch -i sample.csv -o results.csv
```

Optional batch columns include `concentration_mg_ml`, `surface_area_cm2`, `height_cm`, `weight_kg`, `age_years`/`age`, `serum_creatinine_mg_dl`/`creatinine`, and `is_female`/`female`. Missing patient measurements are left missing rather than replaced with assumed values.

## Development

The runtime code uses only the Python standard library.

```bash
python -m pip install pytest build
python -m pip install .
python -m compileall -q chemotherapy_extravasation.py cli.py
python -m pytest -v
python -m build
python cli.py batch -i sample.csv -o out_smoke.csv
```

The CI workflow tests Python 3.9, 3.11, and 3.14, builds the distribution, runs CLI smoke tests, checks installed dependencies, and validates the browser JavaScript and Pages inputs.

## Clinical basis and scope

The current implementation was reviewed against:

- Oncology Nursing Society / American Society of Clinical Oncology, **Extravasation Management: Clinical Practice Guideline**. *JCO Oncology Practice*. DOI: https://doi.org/10.1200/OP-25-00579
- National Cancer Institute, **Common Terminology Criteria for Adverse Events (CTCAE) v5.0**, specifically the “Infusion site extravasation” term.
- eviQ, **Extravasation management – Clinical procedure**, used as a secondary practical cross-check for response workflow and supportive care.

Guidelines, product labeling, institutional policy, available formulations, and drug classifications can change. The source registry is intentionally finite; unknown agents are rejected rather than assigned a default classification or treatment.

## Technology

- Python 3.9+
- Pyodide 314.0.7 (CPython/WebAssembly in the browser)
- HTML, CSS, and vanilla JavaScript
- GitHub Actions and GitHub Pages

## License

MIT. See [LICENSE](LICENSE).
