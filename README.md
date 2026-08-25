# Chemotherapy Extravasation Clinical Decision Support Engine

> **Domain:** Medical Oncology / Chemotherapy Safety & Extravasation Protocol Management  
> **Clinical Guidelines:** ASCO (American Society of Clinical Oncology), ONS (Oncology Nursing Society), and ESMO (European Society for Medical Oncology)

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.9%2B-3776AB.svg?logo=python&logoColor=white)
![Tests](https://img.shields.io/badge/Tests-30%2F30%20Passing-brightgreen.svg)
![Coverage](https://img.shields.io/badge/Coverage-100%25-brightgreen.svg)

---

## Overview

The **Chemotherapy Extravasation Clinical Decision Support Engine** is a deterministic, clinical-grade algorithmic platform designed to guide medical oncologists, oncology nurses, and clinical pharmacists through the acute emergency management and risk assessment of antineoplastic extravasations.

Extravasation of vesicant antineoplastics is an oncologic emergency that can lead to severe tissue necrosis, tendon and nerve damage, compartment syndrome, and permanent loss of extremity function if not managed immediately according to evidence-based pharmacological protocols.

---

## Key Clinical Features

- **Extensive Antineoplastic Registry (25+ Agents):** Categorization into DNA-binding vesicants, non-DNA-binding vesicants, irritants with vesicant potential, irritants, and non-vesicants.
- **Specific Antidote Dosing Engines:**
  - **Dexrazoxane (Totect / Savene):** 3-day systemic IV infusion protocol (1000 mg/m², 1000 mg/m², 500 mg/m²) for anthracycline extravasations within a 6-hour therapeutic window, featuring automatic 50% dose reduction for renal impairment ($\text{CrCl} < 50\text{ mL/min}$).
  - **Hyaluronidase (Amphadase / Vitrase):** Radial clock-face subcutaneous infiltration (150–1500 USP units) for Vinca alkaloids, Taxanes, and Etoposide within a 1-hour window.
  - **Sodium Thiosulfate (STS 10%):** Subcutaneous neutralization protocol (2 mL per 1 mL extravasate, capped at 10 mL) for Cisplatin, Mechlorethamine, and Mitomycin-C.
  - **Dimethyl Sulfoxide (DMSO 99% topical):** Free-radical scavenger protocol for Anthracyclines (when Dexrazoxane is unavailable) or Mitomycin-C.
- **Thermal Protocol Arbitration:**
  - **Dry Cold Compresses:** 15–20 min QID for DNA-binding vesicants, Cisplatin, and Alkylating agents (localizes agent and retards necrosis).
  - **Dry Warm Compresses:** 15–20 min QID for Vinca alkaloids and Etoposide (promotes vasodilation and enzymatic drug clearance). *Cold application is strictly contraindicated for Vinca alkaloids.*
- **CTCAE v5.0 Severity Grading & Surgical Triggers:** Automated classification into Grades 1 to 4 with immediate plastic/orthopedic surgical consult determinations.
- **Multifactorial Risk Scoring (0–100 Scale):** Pre-infusion risk scoring incorporating drug vesicant class, vascular access device type, patient anatomical risk factors, and infusion pump parameters.

---

## Clinical Formulas & Pharmacological Models

### 1. Mosteller Body Surface Area (BSA)
$$\text{BSA } (\text{m}^2) = \sqrt{\frac{\text{Height (cm)} \times \text{Weight (kg)}}{3600}}$$

### 2. Cockcroft-Gault Creatinine Clearance (CrCl)
$$\text{CrCl } (\text{mL/min}) = \frac{(140 - \text{Age}) \times \text{Weight (kg)}}{72 \times \text{Serum Creatinine (mg/dL)}} \times [0.85 \text{ if Female}]$$

### 3. Dexrazoxane Dosing Schedule (Anthracyclines)
$$\text{Day 1 Dose} = \min\left(1000 \times \text{BSA} \times [0.5 \text{ if CrCl} < 50], 2000 \times [0.5 \text{ if CrCl} < 50]\right) \text{ mg}$$
$$\text{Day 2 Dose} = \min\left(1000 \times \text{BSA} \times [0.5 \text{ if CrCl} < 50], 2000 \times [0.5 \text{ if CrCl} < 50]\right) \text{ mg}$$
$$\text{Day 3 Dose} = \min\left(500 \times \text{BSA} \times [0.5 \text{ if CrCl} < 50], 1000 \times [0.5 \text{ if CrCl} < 50]\right) \text{ mg}$$

---

## Command-Line Interface (CLI)

The CLI provides commands for acute event assessment, antidote calculation, risk stratification, staging, and batch CSV processing.

### 1. Comprehensive Extravasation Assessment
```bash
python cli.py assess --drug doxorubicin --catheter peripheral_forearm --volume 8.0 --elapsed 1.5 --height 175 --weight 72 --age 58 --creatinine 1.1 --pain 5
```

### 2. Antidote Regimen Calculator
```bash
# Calculate Dexrazoxane with renal clearance assessment
python cli.py antidote --type dexrazoxane --height 170 --weight 65 --age 62 --creatinine 1.6 --female --elapsed 2.0

# Calculate Hyaluronidase for Vincristine extravasation
python cli.py antidote --type hyaluronidase --volume 4.0 --elapsed 0.5
```

### 3. Pre-Infusion Risk Stratification
```bash
python cli.py risk --drug vincristine --catheter peripheral_hand_wrist --fragile-veins --multiple-attempts
```

### 4. CTCAE v5.0 Severity Grading
```bash
python cli.py stage --pain 8 --ulceration
```

### 5. Batch CSV Processing
```bash
python cli.py batch --input sample.csv --output results.csv
```

### 6. Interactive Clinical Wizard
```bash
python cli.py interactive
```

---

## Python API Usage

```python
from chemotherapy_extravasation import ChemotherapyExtravasationEngine, assess_extravasation_risk

engine = ChemotherapyExtravasationEngine()

# Evaluate acute extravasation event
dossier = engine.evaluate_extravasation_event(
    drug_name="doxorubicin",
    catheter_type="peripheral_forearm",
    estimated_volume_ml=10.0,
    time_elapsed_hours=1.5,
    patient_height_cm=175.0,
    patient_weight_kg=75.0,
    patient_age_years=55,
    serum_creatinine_mg_dl=1.0,
    pain_score_0_to_10=6,
)

print(f"Vesicant Class:    {dossier.vesicant_class}")
print(f"Thermal Protocol:  {dossier.thermal_protocol['protocol']}")
print(f"Antidote Indicated:{dossier.antidote_protocol['antidote_name']}")
print(f"Dosing Summary:    {dossier.antidote_protocol['dose_summary']}")
```

---

## Unit Test Suite

Run the full automated unit test suite:

```bash
python -m unittest test_chemotherapy_extravasation.py
```

Test coverage includes:
- Antineoplastic drug classification and registry invariants
- Thermal compress rules & contraindication verification
- Dexrazoxane BSA calculations, dose caps, and renal CrCl adjustments
- Hyaluronidase volume-based scaling and injection geometries
- Sodium Thiosulfate and DMSO calculations
- CTCAE v5.0 grading and surgical consult triggers
- Multifactorial risk scoring matrix across composite patient and access parameters
- Batch CSV pipeline validation

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
