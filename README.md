# Chemotherapy Extravasation Agent

> **Domain:** Clinical Pharmacology & Precision Pharmacotherapy  
> **Reference Guidelines & Standards:** `CPIC Guidelines & FDA Table of Pharmacogenomic Biomarkers`

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-0.111-009688.svg?logo=fastapi&logoColor=white)
![Audit Trail](https://img.shields.io/badge/Audit-HMAC--SHA256_Tamper--Evident-brightgreen.svg)
![Zero-PHI Guard](https://img.shields.io/badge/Guard-Zero--PHI_Outbound-blue.svg)
![Docker](https://img.shields.io/badge/Docker-Ready-2496ED.svg?logo=docker&logoColor=white)

</div>

---

## 📖 What It Does

**Chemotherapy Extravasation Agent** is an advanced analytical and computational platform implementing Vesicant Extravasation & Specific Antidote (Dexrazoxane) Protocol.

---

## ⚙️ Key Capabilities & Algorithmic Modules

### 🔬 Core Algorithmic & Evaluation Engines

- **`VesicantClass`** — dedicated module for vesicant class evaluation and state verification.
- **`ThermalProtocol`** — dedicated module for thermal protocol evaluation and state verification.
- **`AntidoteType`** — dedicated module for antidote type evaluation and state verification.
- **`CatheterType`** — dedicated module for catheter type evaluation and state verification.
- **`CTCAEGrade`** — dedicated module for c t c a e grade evaluation and state verification.
- **`AntineoplasticDrugInfo`** — dedicated module for antineoplastic drug info evaluation and state verification.

---

## 📐 Mathematical Formulation & Logic

```text
  Mosteller formula for Body Surface Area (BSA):
  Cockcroft-Gault formula for Creatinine Clearance:
  """Calculates specific clinical antidote dosing regimens."""
  bsa = calculate_bsa_mosteller(patient_height_cm, patient_weight_kg)
  antidote_plan = self.antidote_calc.calculate_dmso(25.0, time_elapsed_hours)
```

---

## 💻 CLI Quickstart & Usage

### 1. Guided Interactive Mode
```bash
python cli.py
```

### 2. Direct Parameterized Evaluation
```bash
python cli.py --drug <value> --catheter <value> --volume <value> --elapsed <value>
```

### Parameter Reference
- `--drug`: Specifies input measurement or parameter value.
- `--catheter`: Specifies input measurement or parameter value.
- `--volume`: Specifies input measurement or parameter value.
- `--elapsed`: Specifies input measurement or parameter value.
- `--weight`: Specifies input measurement or parameter value.
- `--height`: Specifies input measurement or parameter value.
- `--age`: Specifies input measurement or parameter value.
- `--creatinine`: Specifies input measurement or parameter value.
- `--type`: Specifies input measurement or parameter value.
- `--bsa`: Specifies input measurement or parameter value.

### Input Data Schema

| Field | Description | Requirement |
|:------|:------------|:------------|
| `case_id` | Parameter / observation metric | Required |
| `patient_id` | Parameter / observation metric | Required |
| `drug` | Parameter / observation metric | Required |
| `catheter` | Parameter / observation metric | Required |
| `volume_ml` | Parameter / observation metric | Required |
| `elapsed_hours` | Parameter / observation metric | Required |
| `height_cm` | Parameter / observation metric | Required |
| `weight_kg` | Parameter / observation metric | Required |

---

## 🛡️ Security & Enterprise Architecture

* **Zero-PHI Outbound Interceptor:** Active AST and regex inspection blocking SSNs, MRNs, phone numbers, and patient identifiers.
* **Tamper-Evident HMAC-SHA256 Audit Trail:** Chained, cryptographically signed logs for every evaluation and state transition.
* **Air-Gapped LLM Reasoning Adapter:** Agnostic integration for local Ollama instances (`llama3`, `mistral`), Claude 3.5 Sonnet, GPT-4o, and deterministic test mocks.
* **Active Learning Bayesian Calibration:** Dynamic tracker updating worker reliability weights and monitoring Brier calibration drift.
* **FastAPI & Prometheus Telemetry:** Exposes OpenAPI 3.1 REST endpoints and operational Prometheus metrics (`/metrics`).

---

## 🧪 Testing & Verification

Run the automated test suite:

```bash
pytest -v
```

Execute high-throughput batch simulation benchmarks:

```bash
python simulator.py --tasks 1000 --concurrency 8
```

---

## 🐳 Container Deployment

```bash
docker build -t chemotherapy-extravasation-agent .
docker run -p 8000:8000 chemotherapy-extravasation-agent
```
