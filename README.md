# Chemotherapy Extravasation Clinical Decision Support Engine

> **Domain:** Medical Oncology, Infusion Nursing & Emergency Antidote Pharmacology  
> **Clinical Guidelines & Standards:** ASCO/ONS Chemotherapy Extravasation Guidelines, ESMO-EONS Clinical Practice Guidelines, CTCAE v5.0, Mosteller BSA Formulation

---

## 📖 Clinical Overview

The **Chemotherapy Extravasation Clinical Decision Support Engine** provides immediate point-of-care risk stratification, CTCAE v5.0 grading, thermal application guidance, and specific pharmacological antidote dosing (Dexrazoxane, Hyaluronidase, Sodium Thiosulfate, DMSO) upon accidental extravasation of antineoplastic agents.

### Vesicant Classification & Antidote Protocols

| Drug Class | Exemplar Agents | Thermal Management | Specific Pharmacological Antidote | Dosing Protocol |
|:---|:---|:---|:---|:---|
| **DNA-Binding Vesicants** | Doxorubicin, Daunorubicin, Epirubicin, Idarubicin | **Dry Cold Pack** (15-20 min q4h) | **Dexrazoxane (Totect / Savene)** | Day 1: $1000\,\text{mg/m}^2$ (max 2000 mg); Day 2: $1000\,\text{mg/m}^2$; Day 3: $500\,\text{mg/m}^2$ (reduce by 50% if $\text{CrCl} < 40\,\text{mL/min}$) |
| **Non-DNA-Binding Vesicants** | Vincristine, Vinblastine, Vinorelbine | **Dry Warm Pack** (15-20 min q4h) | **Hyaluronidase** | 150–1500 units subcutaneous clockwise around extravasation site |
| **Non-DNA-Binding Vesicants (Taxanes)** | Paclitaxel, Docetaxel | **Dry Cold Pack** | **Hyaluronidase** (or cold only) | 150–250 units clockwise around perimeter |
| **Alkylating Agents** | Cisplatin ($\ge 0.5\,\text{mg/mL}$), Mechlorethamine | **Dry Cold Pack** | **Sodium Thiosulfate (1/6 M)** | Infiltrate 2 mL for each 1 mL extravasated |
| **Irritants with Vesicant Potential**| Oxaliplatin, Carboplatin, Mitomycin | **Dry Cold Pack** (Warm for Oxaliplatin) | Drug-specific or symptomatic | Supportive cold/warm compress, site elevation |

---

## 💻 CLI Quickstart & Usage

### 1. Comprehensive Extravasation Assessment
```bash
python cli.py assess --drug doxorubicin --catheter peripheral_forearm --volume 8.0 --elapsed 1.2 --height 172 --weight 68.5 --age 58 --creatinine 1.0 --pain 5
```

### 2. Isolated Antidote Calculation
```bash
python cli.py antidote --drug doxorubicin --height 172 --weight 68.5 --creatinine 1.0 --elapsed 1.5
```

### 3. Interactive Clinical Wizard
```bash
python cli.py interactive
```

### 4. Batch Process Extravasation Incident CSV Log
```bash
python cli.py batch -i sample.csv -o out_results.csv
```

---

## 🧪 Verification & Testing

Execute comprehensive unit tests via pytest:
```bash
python -m pytest -p no:zarr
```
