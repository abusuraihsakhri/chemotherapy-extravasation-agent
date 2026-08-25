#!/usr/bin/env python3
"""
Chemotherapy Extravasation Risk Stratification & Antidote Protocol Engine
========================================================================
Comprehensive clinical decision support module for antineoplastic extravasation
management adhering to ASCO (American Society of Clinical Oncology), ONS (Oncology
Nursing Society), and ESMO (European Society for Medical Oncology) clinical guidelines.

Key Clinical Capabilities:
- Pharmacological drug registry (25+ antineoplastic agents with vesicant classification)
- Specific antidote dosing (Dexrazoxane with renal CrCl adjustments, Hyaluronidase,
  Sodium Thiosulfate, DMSO)
- Thermal compress arbitration (Dry Cold vs Dry Warm protocols)
- CTCAE v5.0 Severity Staging & surgical consultation triggers
- Multifactorial extravasation risk scoring (0-100 scale)
- Comprehensive step-by-step emergency clinical protocols
- Batch CSV analysis and reporting

Stdlib only — no external dependencies.
"""

import csv
import datetime
import math
from dataclasses import dataclass, field, asdict
from enum import Enum
from typing import Dict, Any, List, Optional, Tuple


# ==============================================================================
# ENUMS & CONSTANTS
# ==============================================================================

class VesicantClass(str, Enum):
    DNA_BINDING_VESICANT = "dna_binding_vesicant"
    NON_DNA_BINDING_VESICANT = "non_dna_binding_vesicant"
    IRRITANT_WITH_VESICANT_POTENTIAL = "irritant_with_vesicant_potential"
    IRRITANT = "irritant"
    NON_VESICANT = "non_vesicant"


class ThermalProtocol(str, Enum):
    DRY_COLD = "dry_cold"
    DRY_WARM = "dry_warm"
    NEUTRAL_NONE = "neutral_none"


class AntidoteType(str, Enum):
    DEXRAZOXANE = "dexrazoxane"
    HYALURONIDASE = "hyaluronidase"
    SODIUM_THIOSULFATE = "sodium_thiosulfate"
    DMSO = "dimethyl_sulfoxide"
    NONE = "none"


class CatheterType(str, Enum):
    PERIPHERAL_HAND_WRIST = "peripheral_hand_wrist"
    PERIPHERAL_FOREARM = "peripheral_forearm"
    PERIPHERAL_ANTECUBITAL = "peripheral_antecubital"
    MIDLINE = "midline"
    PICC = "picc"
    TUNNELED_CVC = "tunneled_cvc"
    IMPLANTED_PORT = "implanted_port"


class CTCAEGrade(int, Enum):
    GRADE_1 = 1  # Mild erythema, mild edema, pain 1-3, no skin breakdown
    GRADE_2 = 2  # Moderate erythema, edema, pain 4-6, blister < 1cm, phlebitis
    GRADE_3 = 3  # Severe erythema, ulceration, necrosis, pain 7-10, debridement indicated
    GRADE_4 = 4  # Life-threatening, urgent operative intervention, compartment syndrome


# ==============================================================================
# DRUG KNOWLEDGE BASE
# ==============================================================================

@dataclass(frozen=True)
class AntineoplasticDrugInfo:
    generic_name: str
    brand_names: List[str]
    vesicant_class: VesicantClass
    category: str
    primary_antidote: AntidoteType
    secondary_antidote: AntidoteType
    thermal_protocol: ThermalProtocol
    thermal_rationale: str
    cellular_mechanism: str
    antidote_window_hours: float
    high_concentration_vesicant: bool = False


DRUG_REGISTRY: Dict[str, AntineoplasticDrugInfo] = {
    # Anthracyclines & Antitumor Antibiotics (DNA-Binding Vesicants)
    "doxorubicin": AntineoplasticDrugInfo(
        generic_name="doxorubicin",
        brand_names=["adriamycin", "rubex"],
        vesicant_class=VesicantClass.DNA_BINDING_VESICANT,
        category="Anthracycline",
        primary_antidote=AntidoteType.DEXRAZOXANE,
        secondary_antidote=AntidoteType.DMSO,
        thermal_protocol=ThermalProtocol.DRY_COLD,
        thermal_rationale="Vasoconstriction localizes drug, reduces tissue cellular uptake, and retards metabolic activation.",
        cellular_mechanism="DNA intercalation and topoisomerase II inhibition; retained in tissue creating progressive ulceration.",
        antidote_window_hours=6.0,
    ),
    "daunorubicin": AntineoplasticDrugInfo(
        generic_name="daunorubicin",
        brand_names=["cerubidine", "daunoxome"],
        vesicant_class=VesicantClass.DNA_BINDING_VESICANT,
        category="Anthracycline",
        primary_antidote=AntidoteType.DEXRAZOXANE,
        secondary_antidote=AntidoteType.DMSO,
        thermal_protocol=ThermalProtocol.DRY_COLD,
        thermal_rationale="Vasoconstriction limits cellular exposure and anthracycline complex penetration.",
        cellular_mechanism="DNA intercalation and topoisomerase II inhibition causing severe prolonged necrosis.",
        antidote_window_hours=6.0,
    ),
    "epirubicin": AntineoplasticDrugInfo(
        generic_name="epirubicin",
        brand_names=["ellence", "pharmorubicin"],
        vesicant_class=VesicantClass.DNA_BINDING_VESICANT,
        category="Anthracycline",
        primary_antidote=AntidoteType.DEXRAZOXANE,
        secondary_antidote=AntidoteType.DMSO,
        thermal_protocol=ThermalProtocol.DRY_COLD,
        thermal_rationale="Dry cold application induces localized vasoconstriction and suppresses cellular necrosis.",
        cellular_mechanism="Anthracycline topoisomerase II poisoning with persistent tissue binding.",
        antidote_window_hours=6.0,
    ),
    "idarubicin": AntineoplasticDrugInfo(
        generic_name="idarubicin",
        brand_names=["idamycin"],
        vesicant_class=VesicantClass.DNA_BINDING_VESICANT,
        category="Anthracycline",
        primary_antidote=AntidoteType.DEXRAZOXANE,
        secondary_antidote=AntidoteType.DMSO,
        thermal_protocol=ThermalProtocol.DRY_COLD,
        thermal_rationale="Dry cold induces localized vasoconstriction to restrict tissue distribution.",
        cellular_mechanism="Highly lipophilic DNA-intercalating anthracycline causing deep tissue damage.",
        antidote_window_hours=6.0,
    ),
    "mitomycin": AntineoplasticDrugInfo(
        generic_name="mitomycin",
        brand_names=["mutamycin", "mitomycin-c"],
        vesicant_class=VesicantClass.DNA_BINDING_VESICANT,
        category="Antitumor Antibiotic",
        primary_antidote=AntidoteType.SODIUM_THIOSULFATE,
        secondary_antidote=AntidoteType.DMSO,
        thermal_protocol=ThermalProtocol.DRY_COLD,
        thermal_rationale="Cold compress minimizes tissue perfusion and localizes alkylating metabolite exposure.",
        cellular_mechanism="DNA alkylation and cross-linking leading to persistent, delayed ulceration.",
        antidote_window_hours=2.0,
    ),
    "dactinomycin": AntineoplasticDrugInfo(
        generic_name="dactinomycin",
        brand_names=["cosmegen", "actinomycin-d"],
        vesicant_class=VesicantClass.DNA_BINDING_VESICANT,
        category="Antitumor Antibiotic",
        primary_antidote=AntidoteType.SODIUM_THIOSULFATE,
        secondary_antidote=AntidoteType.NONE,
        thermal_protocol=ThermalProtocol.DRY_COLD,
        thermal_rationale="Dry cold retards local blood flow and reduces cellular extravasation injury.",
        cellular_mechanism="DNA transcription inhibition with severe delayed tissue sloughing.",
        antidote_window_hours=2.0,
    ),

    # Vinca Alkaloids (Non-DNA-Binding Vesicants - WARM COMPRESS MANDATORY)
    "vincristine": AntineoplasticDrugInfo(
        generic_name="vincristine",
        brand_names=["oncovin", "vincasar"],
        vesicant_class=VesicantClass.NON_DNA_BINDING_VESICANT,
        category="Vinca Alkaloid",
        primary_antidote=AntidoteType.HYALURONIDASE,
        secondary_antidote=AntidoteType.NONE,
        thermal_protocol=ThermalProtocol.DRY_WARM,
        thermal_rationale="Vasodilation accelerates systemic absorption and dilution. COLD IS STRICTLY CONTRAINDICATED as it increases local tissue necrosis.",
        cellular_mechanism="Tubulin polymerization inhibition; metabolized rapidly once dispersed.",
        antidote_window_hours=1.0,
    ),
    "vinblastine": AntineoplasticDrugInfo(
        generic_name="vinblastine",
        brand_names=["velban"],
        vesicant_class=VesicantClass.NON_DNA_BINDING_VESICANT,
        category="Vinca Alkaloid",
        primary_antidote=AntidoteType.HYALURONIDASE,
        secondary_antidote=AntidoteType.NONE,
        thermal_protocol=ThermalProtocol.DRY_WARM,
        thermal_rationale="Dry warm compresses promote capillary vasodilation and drug dispersal. Cold application enhances ulceration.",
        cellular_mechanism="Microtubule disassembly arrest in metaphase.",
        antidote_window_hours=1.0,
    ),
    "vinorelbine": AntineoplasticDrugInfo(
        generic_name="vinorelbine",
        brand_names=["navelbine"],
        vesicant_class=VesicantClass.NON_DNA_BINDING_VESICANT,
        category="Vinca Alkaloid",
        primary_antidote=AntidoteType.HYALURONIDASE,
        secondary_antidote=AntidoteType.NONE,
        thermal_protocol=ThermalProtocol.DRY_WARM,
        thermal_rationale="Dry warm heat accelerates clearance through vascular dilation. Avoid cold.",
        cellular_mechanism="Semi-synthetic vinca alkaloid causing localized phlebitis and necrosis if sequestered.",
        antidote_window_hours=1.0,
    ),
    "vindesine": AntineoplasticDrugInfo(
        generic_name="vindesine",
        brand_names=["eldisine"],
        vesicant_class=VesicantClass.NON_DNA_BINDING_VESICANT,
        category="Vinca Alkaloid",
        primary_antidote=AntidoteType.HYALURONIDASE,
        secondary_antidote=AntidoteType.NONE,
        thermal_protocol=ThermalProtocol.DRY_WARM,
        thermal_rationale="Warm compresses facilitate enzymatic interstitial spreading and absorption.",
        cellular_mechanism="Mitotic spindle disruption with non-DNA binding cytotoxicity.",
        antidote_window_hours=1.0,
    ),

    # Taxanes & Trabectedin (Non-DNA-Binding Vesicants / Irritants with Vesicant Properties)
    "paclitaxel": AntineoplasticDrugInfo(
        generic_name="paclitaxel",
        brand_names=["taxol", "abraxane"],
        vesicant_class=VesicantClass.NON_DNA_BINDING_VESICANT,
        category="Taxane",
        primary_antidote=AntidoteType.HYALURONIDASE,
        secondary_antidote=AntidoteType.NONE,
        thermal_protocol=ThermalProtocol.DRY_COLD,
        thermal_rationale="Cold compress reduces local inflammation, Cremophor EL carrier irritation, and edema.",
        cellular_mechanism="Microtubule stabilization preventing depolymerization.",
        antidote_window_hours=1.0,
    ),
    "docetaxel": AntineoplasticDrugInfo(
        generic_name="docetaxel",
        brand_names=["taxotere"],
        vesicant_class=VesicantClass.NON_DNA_BINDING_VESICANT,
        category="Taxane",
        primary_antidote=AntidoteType.HYALURONIDASE,
        secondary_antidote=AntidoteType.NONE,
        thermal_protocol=ThermalProtocol.DRY_COLD,
        thermal_rationale="Cold compresses restrict localized inflammatory mediator release and tissue damage.",
        cellular_mechanism="Microtubule hyperstabilization and mitotic arrest.",
        antidote_window_hours=1.0,
    ),
    "trabectedin": AntineoplasticDrugInfo(
        generic_name="trabectedin",
        brand_names=["yondelis"],
        vesicant_class=VesicantClass.NON_DNA_BINDING_VESICANT,
        category="Alkylating Agent",
        primary_antidote=AntidoteType.DMSO,
        secondary_antidote=AntidoteType.NONE,
        thermal_protocol=ThermalProtocol.DRY_COLD,
        thermal_rationale="Cold compress limits tissue spread of minor groove alkylator.",
        cellular_mechanism="Binds DNA minor groove, bending helix towards major groove.",
        antidote_window_hours=2.0,
    ),

    # Platinum Agents & Alkylating Irritants with Vesicant Potential
    "cisplatin": AntineoplasticDrugInfo(
        generic_name="cisplatin",
        brand_names=["platinol"],
        vesicant_class=VesicantClass.IRRITANT_WITH_VESICANT_POTENTIAL,
        category="Platinum Coordination Complex",
        primary_antidote=AntidoteType.SODIUM_THIOSULFATE,
        secondary_antidote=AntidoteType.NONE,
        thermal_protocol=ThermalProtocol.DRY_COLD,
        thermal_rationale="Cold compress localizes platinum ion diffusion. Apply STS if > 20 mL of >= 0.5 mg/mL extravasates.",
        cellular_mechanism="Intrastrand and interstrand DNA crosslinking.",
        antidote_window_hours=2.0,
        high_concentration_vesicant=True,
    ),
    "carboplatin": AntineoplasticDrugInfo(
        generic_name="carboplatin",
        brand_names=["paraplatin"],
        vesicant_class=VesicantClass.IRRITANT,
        category="Platinum Coordination Complex",
        primary_antidote=AntidoteType.SODIUM_THIOSULFATE,
        secondary_antidote=AntidoteType.NONE,
        thermal_protocol=ThermalProtocol.DRY_COLD,
        thermal_rationale="Cold compress limits pain, swelling, and chemical phlebitis.",
        cellular_mechanism="DNA adduct formation with lower local vesicant reactivity than cisplatin.",
        antidote_window_hours=2.0,
    ),
    "oxaliplatin": AntineoplasticDrugInfo(
        generic_name="oxaliplatin",
        brand_names=["eloxatin"],
        vesicant_class=VesicantClass.IRRITANT_WITH_VESICANT_POTENTIAL,
        category="Platinum Coordination Complex",
        primary_antidote=AntidoteType.SODIUM_THIOSULFATE,
        secondary_antidote=AntidoteType.NONE,
        thermal_protocol=ThermalProtocol.DRY_WARM,
        thermal_rationale="WARM or room temp compress advised because cold exposure can trigger acute neurosensory dysesthesia / pharyngolaryngeal spasm.",
        cellular_mechanism="Platinum-DNA crosslinks; neurotoxicity hypersensitive to cold.",
        antidote_window_hours=2.0,
    ),
    "dacarbazine": AntineoplasticDrugInfo(
        generic_name="dacarbazine",
        brand_names=["dtic-dome"],
        vesicant_class=VesicantClass.IRRITANT_WITH_VESICANT_POTENTIAL,
        category="Alkylating Agent",
        primary_antidote=AntidoteType.SODIUM_THIOSULFATE,
        secondary_antidote=AntidoteType.NONE,
        thermal_protocol=ThermalProtocol.DRY_COLD,
        thermal_rationale="Cold compress minimizes tissue exposure to photoreactive alkylating degradation products.",
        cellular_mechanism="DNA methylation through diazomethane active intermediate.",
        antidote_window_hours=2.0,
    ),
    "carmustine": AntineoplasticDrugInfo(
        generic_name="carmustine",
        brand_names=["bicnu"],
        vesicant_class=VesicantClass.IRRITANT_WITH_VESICANT_POTENTIAL,
        category="Nitrosourea",
        primary_antidote=AntidoteType.SODIUM_THIOSULFATE,
        secondary_antidote=AntidoteType.NONE,
        thermal_protocol=ThermalProtocol.DRY_COLD,
        thermal_rationale="Cold compresses restrict nitrosourea lipophilic tissue penetration.",
        cellular_mechanism="DNA and RNA carbamoylation and crosslinking.",
        antidote_window_hours=2.0,
    ),
    "ifosfamide": AntineoplasticDrugInfo(
        generic_name="ifosfamide",
        brand_names=["ifex"],
        vesicant_class=VesicantClass.IRRITANT,
        category="Alkylating Agent",
        primary_antidote=AntidoteType.NONE,
        secondary_antidote=AntidoteType.NONE,
        thermal_protocol=ThermalProtocol.DRY_COLD,
        thermal_rationale="Cold compress relieves localized burning and erythema.",
        cellular_mechanism="Alkylating oxazaphosphorine requiring hepatic activation.",
        antidote_window_hours=2.0,
    ),
    "etoposide": AntineoplasticDrugInfo(
        generic_name="etoposide",
        brand_names=["vepesid", "toposar"],
        vesicant_class=VesicantClass.IRRITANT_WITH_VESICANT_POTENTIAL,
        category="Topoisomerase II Inhibitor",
        primary_antidote=AntidoteType.HYALURONIDASE,
        secondary_antidote=AntidoteType.NONE,
        thermal_protocol=ThermalProtocol.DRY_WARM,
        thermal_rationale="Dry warm compresses promote absorption and dispersion of concentrated solution.",
        cellular_mechanism="Topoisomerase II complex stabilization preventing DNA religation.",
        antidote_window_hours=1.0,
        high_concentration_vesicant=True,
    ),

    # Common Antimetabolites & Alkylating Irritants / Non-Vesicants
    "fluorouracil": AntineoplasticDrugInfo(
        generic_name="fluorouracil",
        brand_names=["5-fu", "adrucil"],
        vesicant_class=VesicantClass.IRRITANT,
        category="Antimetabolite (Pyrimidine)",
        primary_antidote=AntidoteType.NONE,
        secondary_antidote=AntidoteType.NONE,
        thermal_protocol=ThermalProtocol.DRY_COLD,
        thermal_rationale="Cold compress alleviates acute hyperpigmentation, phlebitis, and local burning.",
        cellular_mechanism="Thymidylate synthase inhibition and RNA misincorporation.",
        antidote_window_hours=0.0,
    ),
    "gemcitabine": AntineoplasticDrugInfo(
        generic_name="gemcitabine",
        brand_names=["gemzar"],
        vesicant_class=VesicantClass.IRRITANT,
        category="Antimetabolite (Pyrimidine)",
        primary_antidote=AntidoteType.NONE,
        secondary_antidote=AntidoteType.NONE,
        thermal_protocol=ThermalProtocol.DRY_COLD,
        thermal_rationale="Cold compresses mitigate erythema and inflammatory induration.",
        cellular_mechanism="DNA polymerase inhibition and ribonucleotide reductase blockade.",
        antidote_window_hours=0.0,
    ),
    "cyclophosphamide": AntineoplasticDrugInfo(
        generic_name="cyclophosphamide",
        brand_names=["cytoxan", "neosar"],
        vesicant_class=VesicantClass.IRRITANT,
        category="Alkylating Agent (Nitrogen Mustard)",
        primary_antidote=AntidoteType.NONE,
        secondary_antidote=AntidoteType.NONE,
        thermal_protocol=ThermalProtocol.DRY_COLD,
        thermal_rationale="Cold compress reduces local discomfort. Requires liver activation so non-toxic locally.",
        cellular_mechanism="Hepatic prodrug metabolized to phosphoramide mustard.",
        antidote_window_hours=0.0,
    ),
    "methotrexate": AntineoplasticDrugInfo(
        generic_name="methotrexate",
        brand_names=["trexall", "otrexup"],
        vesicant_class=VesicantClass.IRRITANT,
        category="Antimetabolite (Folate Antagonist)",
        primary_antidote=AntidoteType.NONE,
        secondary_antidote=AntidoteType.NONE,
        thermal_protocol=ThermalProtocol.DRY_COLD,
        thermal_rationale="Cold compress reduces localized edema and soreness.",
        cellular_mechanism="Dihydrofolate reductase (DHFR) competitive inhibition.",
        antidote_window_hours=0.0,
    ),
    "bleomycin": AntineoplasticDrugInfo(
        generic_name="bleomycin",
        brand_names=["blenoxane"],
        vesicant_class=VesicantClass.NON_VESICANT,
        category="Antitumor Antibiotic (Glycopeptide)",
        primary_antidote=AntidoteType.NONE,
        secondary_antidote=AntidoteType.NONE,
        thermal_protocol=ThermalProtocol.DRY_COLD,
        thermal_rationale="Cold compress for symptom relief.",
        cellular_mechanism="DNA scission via iron-dependent oxygen radical complex.",
        antidote_window_hours=0.0,
    ),
}


# ==============================================================================
# CLINICAL CALCULATION MODELS
# ==============================================================================

def calculate_bsa_mosteller(height_cm: float, weight_kg: float) -> float:
    """
    Mosteller formula for Body Surface Area (BSA):
    BSA (m²) = sqrt((height_cm * weight_kg) / 3600)
    """
    if height_cm <= 0 or weight_kg <= 0:
        raise ValueError("Height and weight must be positive numbers.")
    bsa = math.sqrt((height_cm * weight_kg) / 3600.0)
    return round(bsa, 2)


def calculate_crcl_cockcroft_gault(age_years: int, weight_kg: float, serum_creatinine_mg_dl: float, is_female: bool) -> float:
    """
    Cockcroft-Gault formula for Creatinine Clearance:
    CrCl = ((140 - age) * weight_kg) / (72 * Scr) [* 0.85 if female]
    """
    if age_years < 18 or age_years > 120:
        raise ValueError("Age must be between 18 and 120 years.")
    if weight_kg <= 0:
        raise ValueError("Weight must be positive.")
    if serum_creatinine_mg_dl <= 0:
        raise ValueError("Serum creatinine must be positive.")

    crcl = ((140.0 - age_years) * weight_kg) / (72.0 * serum_creatinine_mg_dl)
    if is_female:
        crcl *= 0.85
    return round(crcl, 1)


# ==============================================================================
# ANTIDOTE DOSING CALCULATOR
# ==============================================================================

@dataclass
class AntidotePlan:
    antidote_name: str
    is_indicated: bool
    urgency_window_hours: float
    time_elapsed_hours: float
    is_within_window: bool
    dose_summary: str
    schedule: List[Dict[str, Any]]
    administration_instructions: List[str]
    contraindications_and_warnings: List[str]
    renal_adjustment_applied: bool = False


class AntidoteCalculator:
    """Calculates specific clinical antidote dosing regimens."""

    @staticmethod
    def calculate_dexrazoxane(
        bsa_m2: float,
        crcl_ml_min: Optional[float] = None,
        time_elapsed_hours: float = 0.0,
    ) -> AntidotePlan:
        """
        Dexrazoxane (Totect/Savene) 3-day IV infusion protocol for anthracycline extravasations.
        Dosing:
        - Day 1: 1000 mg/m² (max 2000 mg) within 6 hours.
        - Day 2: 1000 mg/m² (max 2000 mg) at 24 hours ± 3h.
        - Day 3: 500 mg/m² (max 1000 mg) at 48 hours ± 3h.
        Renal adjustment: CrCl < 50 mL/min -> Reduce dose by 50%.
        """
        if bsa_m2 <= 0:
            raise ValueError("BSA must be positive.")

        renal_adj = crcl_ml_min is not None and crcl_ml_min < 50.0
        dose_factor = 0.5 if renal_adj else 1.0

        d1_rate = 1000.0 * dose_factor
        d2_rate = 1000.0 * dose_factor
        d3_rate = 500.0 * dose_factor

        d1_max = 2000.0 * dose_factor
        d2_max = 2000.0 * dose_factor
        d3_max = 1000.0 * dose_factor

        d1_dose = min(round(d1_rate * bsa_m2), d1_max)
        d2_dose = min(round(d2_rate * bsa_m2), d2_max)
        d3_dose = min(round(d3_rate * bsa_m2), d3_max)

        schedule = [
            {
                "day": 1,
                "timing": "Immediately (within 6 hours of extravasation)",
                "dose_mg": d1_dose,
                "dose_per_m2": d1_rate,
                "infusion_duration_minutes": 60,
                "diluent": "500 mL 0.9% Sodium Chloride or D5W",
            },
            {
                "day": 2,
                "timing": "24 hours (±3 hours) after Day 1 dose",
                "dose_mg": d2_dose,
                "dose_per_m2": d2_rate,
                "infusion_duration_minutes": 60,
                "diluent": "500 mL 0.9% Sodium Chloride or D5W",
            },
            {
                "day": 3,
                "timing": "48 hours (±3 hours) after Day 1 dose",
                "dose_mg": d3_dose,
                "dose_per_m2": d3_rate,
                "infusion_duration_minutes": 60,
                "diluent": "500 mL 0.9% Sodium Chloride or D5W",
            },
        ]

        instructions = [
            "Initiate Day 1 infusion as soon as possible and NO LATER than 6 hours post-extravasation.",
            "Infuse intravenously over 1 to 2 hours in a DIFFERENT extremity/vein from the extravasation site.",
            "Do NOT use inline filters during infusion.",
            "Remove cold compresses at least 15 minutes before and during Dexrazoxane infusion to ensure adequate vascular delivery.",
            "Monitor complete blood count (CBC) and liver function tests (LFTs) due to potential transient myelosuppression and transaminitis.",
        ]

        warnings = []
        if renal_adj:
            warnings.append(f"Renal impairment detected (CrCl {crcl_ml_min:.1f} mL/min < 50 mL/min): Dexrazoxane doses reduced by 50%.")
        if time_elapsed_hours > 6.0:
            warnings.append(f"CRITICAL WARNING: Time elapsed ({time_elapsed_hours:.1f}h) exceeds the validated 6-hour therapeutic window. Efficacy may be significantly degraded.")

        return AntidotePlan(
            antidote_name="Dexrazoxane (Totect / Savene)",
            is_indicated=True,
            urgency_window_hours=6.0,
            time_elapsed_hours=time_elapsed_hours,
            is_within_window=time_elapsed_hours <= 6.0,
            dose_summary=f"Day 1: {d1_dose} mg, Day 2: {d2_dose} mg, Day 3: {d3_dose} mg IV (Total: {d1_dose + d2_dose + d3_dose} mg)",
            schedule=schedule,
            administration_instructions=instructions,
            contraindications_and_warnings=warnings,
            renal_adjustment_applied=renal_adj,
        )

    @staticmethod
    def calculate_hyaluronidase(
        estimated_volume_ml: float = 5.0,
        time_elapsed_hours: float = 0.0,
    ) -> AntidotePlan:
        """
        Hyaluronidase subcutaneous radial injection protocol for Vinca alkaloids, Taxanes, and Etoposide.
        Dosing: 150 to 300 USP units (1-2 mL of 150 units/mL solution).
        Technique: Clock-face subcutaneous infiltration around circumference.
        """
        if estimated_volume_ml <= 0:
            raise ValueError("Extravasation volume must be positive.")

        if estimated_volume_ml <= 2.0:
            total_units = 150
            num_injections = 4
        elif estimated_volume_ml <= 10.0:
            total_units = 300
            num_injections = 6
        else:
            # For extensive extravasation (>10 mL)
            total_units = min(int(math.ceil(estimated_volume_ml / 5.0) * 150), 1500)
            num_injections = 8

        units_per_injection = round(total_units / num_injections, 1)
        ml_per_injection = round((total_units / 150.0) / num_injections, 2)

        schedule = [
            {
                "step": "Subcutaneous Infiltration",
                "total_units": total_units,
                "num_sites": num_injections,
                "units_per_site": units_per_injection,
                "volume_per_site_ml": ml_per_injection,
                "timing": "Immediately within 1 hour of extravasation",
            }
        ]

        instructions = [
            "Administer within 1 hour of extravasation (most effective if given within 30 minutes).",
            f"Use a 25-gauge to 27-gauge needle; inject {units_per_injection} units ({ml_per_injection} mL) subcutaneously at {num_injections} sites around the leading edge of the extravasation in a clockwise ring.",
            "Change the needle between each injection site to avoid tracking extravasated drug into undamaged tissue.",
            "Also inject 150 units through the retained catheter prior to catheter removal if possible.",
            "Apply DRY WARM compresses for 20 minutes QID for 24-48 hours. DO NOT APPLY COLD.",
        ]

        warnings = []
        if time_elapsed_hours > 1.0:
            warnings.append(f"WARNING: Time elapsed ({time_elapsed_hours:.1f}h) exceeds optimal 1-hour window for enzymatic spreading.")

        return AntidotePlan(
            antidote_name="Hyaluronidase (Amphadase / Vitrase / Hylenex)",
            is_indicated=True,
            urgency_window_hours=1.0,
            time_elapsed_hours=time_elapsed_hours,
            is_within_window=time_elapsed_hours <= 1.0,
            dose_summary=f"{total_units} USP Units SC divided into {num_injections} clockwise radial subcutaneous injections ({units_per_injection} units/site)",
            schedule=schedule,
            administration_instructions=instructions,
            contraindications_and_warnings=warnings,
        )

    @staticmethod
    def calculate_sodium_thiosulfate(
        estimated_volume_ml: float = 5.0,
        time_elapsed_hours: float = 0.0,
    ) -> AntidotePlan:
        """
        Sodium Thiosulfate (STS) neutralization protocol for Cisplatin, Mechlorethamine, Mitomycin-C.
        Dosing: 2 mL of 1/6 M (approx 4.17%) to 10% solution per 1 mL of estimated extravasate (max 10 mL total).
        """
        if estimated_volume_ml <= 0:
            raise ValueError("Extravasation volume must be positive.")

        sts_volume_ml = min(round(2.0 * estimated_volume_ml, 1), 10.0)
        num_injections = 5
        ml_per_site = round(sts_volume_ml / num_injections, 2)

        schedule = [
            {
                "step": "Subcutaneous Neutralization Infiltration",
                "total_volume_ml": sts_volume_ml,
                "solution_concentration": "1/6 M (approx 4.17%) or 10% solution",
                "num_sites": num_injections,
                "volume_per_site_ml": ml_per_site,
            }
        ]

        instructions = [
            "Preparation: Dilute 4 mL of 25% Sodium Thiosulfate with 6 mL of Sterile Water for Injection to yield 10 mL of 10% solution (or use commercially prepared 1/6 M solution).",
            f"Inject {ml_per_site} mL subcutaneously at {num_injections} sites surrounding the extravasation site using a 25-27G needle.",
            "If catheter is still in place, instill 1-2 mL of STS through the catheter hub before removal.",
            "Apply dry cold compresses for 20 minutes QID for 24-48 hours.",
        ]

        warnings = []
        if time_elapsed_hours > 2.0:
            warnings.append(f"WARNING: Time elapsed ({time_elapsed_hours:.1f}h) exceeds optimal 2-hour neutralization window.")

        return AntidotePlan(
            antidote_name="Sodium Thiosulfate (1/6 M or 10%)",
            is_indicated=True,
            urgency_window_hours=2.0,
            time_elapsed_hours=time_elapsed_hours,
            is_within_window=time_elapsed_hours <= 2.0,
            dose_summary=f"{sts_volume_ml} mL of 10% solution SC divided across {num_injections} subcutaneous sites ({ml_per_site} mL/site)",
            schedule=schedule,
            administration_instructions=instructions,
            contraindications_and_warnings=warnings,
        )

    @staticmethod
    def calculate_dmso(
        surface_area_cm2: float = 25.0,
        time_elapsed_hours: float = 0.0,
    ) -> AntidotePlan:
        """
        Dimethyl Sulfoxide (DMSO 99% topical) free-radical scavenger protocol.
        Used for Anthracyclines (when Dexrazoxane is unavailable) or Mitomycin-C.
        """
        drops_required = max(int(round((surface_area_cm2 / 10.0) * 4.0)), 4)

        schedule = [
            {
                "frequency": "Every 8 hours for 7 to 14 consecutive days",
                "dosage_drops": drops_required,
                "area_coverage_cm2": surface_area_cm2,
            }
        ]

        instructions = [
            f"Apply {drops_required} drops of 99% Dimethyl Sulfoxide (DMSO) topically to twice the affected surface area ({surface_area_cm2 * 2:.0f} cm²).",
            "Allow the solution to air dry completely without occlusive dressings or bandages.",
            "Repeat application every 8 hours for a minimum of 7 days, extending to 14 days if symptoms persist.",
            "Avoid combining topical DMSO with Dexrazoxane as DMSO may decrease the systemic efficacy of Dexrazoxane.",
        ]

        warnings = [
            "Patients will experience a characteristic garlic-like breath odor and taste during topical DMSO therapy.",
            "Do NOT apply occlusive dressings over DMSO-treated areas.",
        ]

        return AntidotePlan(
            antidote_name="Dimethyl Sulfoxide (DMSO 99% Topical)",
            is_indicated=True,
            urgency_window_hours=2.0,
            time_elapsed_hours=time_elapsed_hours,
            is_within_window=time_elapsed_hours <= 24.0,
            dose_summary=f"{drops_required} drops of 99% DMSO applied topically TID for 7-14 days",
            schedule=schedule,
            administration_instructions=instructions,
            contraindications_and_warnings=warnings,
        )


# ==============================================================================
# MULTIFACTORIAL RISK SCORER & CTCAE STAGER
# ==============================================================================

@dataclass
class ExtravasationRiskAssessment:
    composite_risk_score: float  # 0 to 100
    risk_tier: str  # LOW, MODERATE, HIGH, CRITICAL
    drug_name: str
    vesicant_class: str
    catheter_type: str
    score_breakdown: Dict[str, float]
    clinical_flags: List[str]
    prevention_recommendations: List[str]


def assess_extravasation_risk(
    drug_name: str,
    catheter_type: str,
    multiple_venipuncture_attempts: bool = False,
    is_elderly_or_fragile_veins: bool = False,
    has_sensory_neuropathy: bool = False,
    prior_radiation_to_limb: bool = False,
    has_lymphedema: bool = False,
    is_agitated_or_confused: bool = False,
    is_infusion_pump_high_pressure: bool = False,
    prolonged_infusion_gt_4h: bool = False,
) -> ExtravasationRiskAssessment:
    """
    Multifactorial Extravasation Risk Scoring Model (0 - 100 scale).
    """
    drug_key = drug_name.strip().lower()
    drug_info = DRUG_REGISTRY.get(drug_key)

    # 1. Drug Vesicant Potency (0 - 35 points)
    if drug_info:
        if drug_info.vesicant_class == VesicantClass.DNA_BINDING_VESICANT:
            drug_pts = 35.0
        elif drug_info.vesicant_class == VesicantClass.NON_DNA_BINDING_VESICANT:
            drug_pts = 28.0
        elif drug_info.vesicant_class == VesicantClass.IRRITANT_WITH_VESICANT_POTENTIAL:
            drug_pts = 20.0
        elif drug_info.vesicant_class == VesicantClass.IRRITANT:
            drug_pts = 10.0
        else:
            drug_pts = 2.0
    else:
        drug_pts = 12.0  # Unknown drug default

    # 2. Catheter & Site Risk (0 - 25 points)
    catheter_map = {
        CatheterType.PERIPHERAL_HAND_WRIST.value: 25.0,
        CatheterType.PERIPHERAL_ANTECUBITAL.value: 20.0,
        CatheterType.PERIPHERAL_FOREARM.value: 14.0,
        CatheterType.MIDLINE.value: 10.0,
        CatheterType.PICC.value: 5.0,
        CatheterType.TUNNELED_CVC.value: 3.0,
        CatheterType.IMPLANTED_PORT.value: 2.0,
    }
    catheter_pts = catheter_map.get(catheter_type.lower(), 15.0)

    # 3. Patient Anatomical & Clinical Factors (0 - 25 points)
    patient_pts = 0.0
    flags = []
    if is_elderly_or_fragile_veins:
        patient_pts += 6.0
        flags.append("Fragile or sclerotic peripheral vasculature.")
    if has_sensory_neuropathy:
        patient_pts += 7.0
        flags.append("Sensory neuropathy impairing early extravasation pain detection.")
    if prior_radiation_to_limb:
        patient_pts += 6.0
        flags.append("Prior radiation field with impaired lymphatic clearance and fibrotic tissue.")
    if has_lymphedema:
        patient_pts += 6.0
        flags.append("Ipsilateral lymphedema compromising lymphatic drainage.")
    if is_agitated_or_confused:
        patient_pts += 5.0
        flags.append("Patient agitation/confusion increasing risk of accidental catheter dislodgement.")
    patient_pts = min(patient_pts, 25.0)

    # 4. Infusion / Technique Factors (0 - 15 points)
    technique_pts = 0.0
    if multiple_venipuncture_attempts:
        technique_pts += 8.0
        flags.append("Multiple venipuncture attempts distal to infusion site.")
    if is_infusion_pump_high_pressure:
        technique_pts += 5.0
        flags.append("High-pressure electronic infusion pump without extravasation sensor.")
    if prolonged_infusion_gt_4h:
        technique_pts += 4.0
        flags.append("Prolonged infusion (>4 hours) increasing vein wall microvascular wear.")
    technique_pts = min(technique_pts, 15.0)

    total_score = round(min(drug_pts + catheter_pts + patient_pts + technique_pts, 100.0), 1)

    if total_score >= 70.0:
        tier = "CRITICAL"
    elif total_score >= 45.0:
        tier = "HIGH"
    elif total_score >= 25.0:
        tier = "MODERATE"
    else:
        tier = "LOW"

    # Prevention recommendations
    recs = []
    if tier in ("CRITICAL", "HIGH"):
        recs.append("MANDATORY: Continuous nursing bedside visualization during vesicant administration.")
        recs.append("Perform free-flowing blood return check every 2 to 3 mL during bolus injection or every 10 minutes during short infusion.")
        if catheter_pts >= 20.0:
            recs.append("STRONGLY RECOMMEND: Secure central venous access (PICC or Port) prior to next cycle.")
    if drug_info and drug_info.vesicant_class in (VesicantClass.DNA_BINDING_VESICANT, VesicantClass.NON_DNA_BINDING_VESICANT):
        recs.append(f"Pre-stage emergency extravasation kit containing {drug_info.primary_antidote.value.upper()} and appropriate thermal pack at bedside.")
    if has_sensory_neuropathy or is_agitated_or_confused:
        recs.append("Implement Q15-minute palpation and skin inspection protocol due to reduced patient symptom self-reporting.")

    return ExtravasationRiskAssessment(
        composite_risk_score=total_score,
        risk_tier=tier,
        drug_name=drug_name,
        vesicant_class=drug_info.vesicant_class.value if drug_info else "unknown",
        catheter_type=catheter_type,
        score_breakdown={
            "drug_potency_score": drug_pts,
            "catheter_site_score": catheter_pts,
            "patient_factors_score": patient_pts,
            "technique_factors_score": technique_pts,
        },
        clinical_flags=flags,
        prevention_recommendations=recs,
    )


def grade_ctcae_severity(
    pain_score_0_to_10: int,
    erythema_present: bool = True,
    edema_present: bool = True,
    blistering_present: bool = False,
    blister_size_cm: float = 0.0,
    ulceration_or_necrosis_present: bool = False,
    tissue_sloughing_or_eschar: bool = False,
    compartment_syndrome_signs: bool = False,
    loss_of_extremity_function: bool = False,
) -> Dict[str, Any]:
    """
    CTCAE v5.0 Extravasation Severity Grading & Surgical Consult Classifier.
    """
    if compartment_syndrome_signs or loss_of_extremity_function:
        grade = CTCAEGrade.GRADE_4
        desc = "Grade 4 (Life-threatening / Disabling): Urgent operative intervention indicated; compartment syndrome or neurovascular compromise."
        surgical_consult = "EMERGENT: Immediate orthopedic / plastic surgery evaluation for emergency fasciotomy / operative exploration."
    elif ulceration_or_necrosis_present or tissue_sloughing_or_eschar or blister_size_cm >= 1.0 or pain_score_0_to_10 >= 7:
        grade = CTCAEGrade.GRADE_3
        desc = "Grade 3 (Severe): Severe pain, ulceration, gross tissue necrosis, or extensive bullae; surgical debridement indicated."
        surgical_consult = "URGENT: Plastic surgery or wound care consultation within 24 hours for evaluation of debridement."
    elif blistering_present or pain_score_0_to_10 >= 4 or (erythema_present and edema_present):
        grade = CTCAEGrade.GRADE_2
        desc = "Grade 2 (Moderate): Moderate erythema and edema with pain 4-6; localized small blisters (< 1 cm); phlebitis."
        surgical_consult = "ADVISORY: Wound care follow-up; request surgical consult if no clinical improvement within 48-72 hours."
    else:
        grade = CTCAEGrade.GRADE_1
        desc = "Grade 1 (Mild): Mild erythema and localized edema with minimal discomfort (pain 1-3); no skin breakdown."
        surgical_consult = "ROUTINE: Outpatient oncology nursing monitoring; surgery consult not immediately required."

    return {
        "ctcae_grade": grade.value,
        "description": desc,
        "surgical_consultation_status": surgical_consult,
        "pain_score": pain_score_0_to_10,
        "requires_urgent_surgical_review": grade.value >= 3,
    }


# ==============================================================================
# EMERGENCY MANAGEMENT PROTOCOL GENERATOR
# ==============================================================================

@dataclass
class ExtravasationEmergencyDossier:
    event_id: str
    patient_id: str
    drug_name: str
    vesicant_class: str
    catheter_type: str
    ctcae_severity: Dict[str, Any]
    thermal_protocol: Dict[str, Any]
    antidote_protocol: Optional[Dict[str, Any]]
    ordered_action_checklist: List[Dict[str, Any]]
    monitoring_and_followup_schedule: List[str]
    documentation_requirements: List[str]
    timestamp_utc: str = ""

    def __post_init__(self):
        if not self.timestamp_utc:
            self.timestamp_utc = datetime.datetime.now(datetime.timezone.utc).isoformat()


class ChemotherapyExtravasationEngine:
    """Master Clinical Orchestrator for Chemotherapy Extravasation Management."""

    def __init__(self):
        self.antidote_calc = AntidoteCalculator()

    def evaluate_extravasation_event(
        self,
        drug_name: str,
        catheter_type: str,
        estimated_volume_ml: float,
        time_elapsed_hours: float,
        patient_height_cm: Optional[float] = None,
        patient_weight_kg: Optional[float] = None,
        patient_age_years: Optional[int] = None,
        serum_creatinine_mg_dl: Optional[float] = None,
        is_female: bool = False,
        pain_score_0_to_10: int = 4,
        blistering_present: bool = False,
        ulceration_or_necrosis_present: bool = False,
        compartment_syndrome_signs: bool = False,
        patient_id: str = "PATIENT-ONC-001",
    ) -> ExtravasationEmergencyDossier:
        drug_key = drug_name.strip().lower()
        drug_info = DRUG_REGISTRY.get(drug_key)

        if not drug_info:
            vesicant_class = VesicantClass.IRRITANT.value
            thermal_mode = ThermalProtocol.DRY_COLD.value
            thermal_rat = "Standard default cold application for unclassified cytotoxic agent."
            primary_antidote = AntidoteType.NONE
            antidote_window = 1.0
        else:
            vesicant_class = drug_info.vesicant_class.value
            thermal_mode = drug_info.thermal_protocol.value
            thermal_rat = drug_info.thermal_rationale
            primary_antidote = drug_info.primary_antidote
            antidote_window = drug_info.antidote_window_hours

        # 1. CTCAE Severity
        ctcae = grade_ctcae_severity(
            pain_score_0_to_10=pain_score_0_to_10,
            blistering_present=blistering_present,
            ulceration_or_necrosis_present=ulceration_or_necrosis_present,
            compartment_syndrome_signs=compartment_syndrome_signs,
        )

        # 2. Antidote Calculation
        antidote_plan = None
        if primary_antidote == AntidoteType.DEXRAZOXANE:
            bsa = 1.73  # Standard adult default
            if patient_height_cm and patient_weight_kg:
                bsa = calculate_bsa_mosteller(patient_height_cm, patient_weight_kg)
            crcl = None
            if patient_age_years and patient_weight_kg and serum_creatinine_mg_dl:
                crcl = calculate_crcl_cockcroft_gault(patient_age_years, patient_weight_kg, serum_creatinine_mg_dl, is_female)
            antidote_plan = self.antidote_calc.calculate_dexrazoxane(bsa, crcl, time_elapsed_hours)

        elif primary_antidote == AntidoteType.HYALURONIDASE:
            antidote_plan = self.antidote_calc.calculate_hyaluronidase(estimated_volume_ml, time_elapsed_hours)

        elif primary_antidote == AntidoteType.SODIUM_THIOSULFATE:
            antidote_plan = self.antidote_calc.calculate_sodium_thiosulfate(estimated_volume_ml, time_elapsed_hours)

        elif primary_antidote == AntidoteType.DMSO:
            antidote_plan = self.antidote_calc.calculate_dmso(25.0, time_elapsed_hours)

        # 3. Thermal Protocol Instructions
        if thermal_mode == ThermalProtocol.DRY_COLD.value:
            thermal_dict = {
                "protocol": "DRY COLD COMPRESS",
                "frequency": "Apply for 15 to 20 minutes every 6 hours (4 times daily) for 24 to 48 hours.",
                "rationale": thermal_rat,
                "warning": "Ensure compress is DRY (place ice pack in sealed plastic bag wrapped in dry towel). DO NOT APPLY MOIST COLD (prevents tissue maceration). If Dexrazoxane is administered, remove cold compress 15 minutes before and during infusion.",
            }
        elif thermal_mode == ThermalProtocol.DRY_WARM.value:
            thermal_dict = {
                "protocol": "DRY WARM COMPRESS",
                "frequency": "Apply for 15 to 20 minutes every 6 hours (4 times daily) for 24 to 48 hours.",
                "rationale": thermal_rat,
                "warning": "CRITICAL: Dry warm only. STRICTLY AVOID COLD COMPRESSES (cold enhances Vinca alkaloid-induced tissue ulceration). Use heating pad or warm dry pack wrapped in cloth (temperature <= 40°C / 104°F).",
            }
        else:
            thermal_dict = {
                "protocol": "ROOM TEMPERATURE / NO THERMAL INTERVENTION",
                "frequency": "Maintain extremity at ambient temperature.",
                "rationale": thermal_rat,
                "warning": "Do not apply extreme hot or cold packs.",
            }

        # 4. Ordered Step-by-Step Action Checklist
        actions = [
            {
                "step_number": 1,
                "priority": "IMMEDIATE (STAT)",
                "action": "STOP INFUSION IMMEDIATELY",
                "details": "Immediately halt the antineoplastic infusion upon first suspicion or report of pain, burning, swelling, or loss of blood return.",
            },
            {
                "step_number": 2,
                "priority": "IMMEDIATE (STAT)",
                "action": "LEAVE CATHETER IN SITU; DO NOT FLUSH",
                "details": "Disconnect IV administration set from catheter hub. DO NOT FLUSH the catheter with saline or heparin under any circumstances.",
            },
            {
                "step_number": 3,
                "priority": "IMMEDIATE (STAT)",
                "action": "ASPIRATE EXTRAVASATED DRUG & BLOOD",
                "details": "Attach a sterile 3 mL or 5 mL syringe to the catheter hub and gently aspirate 3 to 5 mL of blood and extravasated fluid.",
            },
        ]

        if primary_antidote in (AntidoteType.SODIUM_THIOSULFATE, AntidoteType.HYALURONIDASE):
            actions.append({
                "step_number": 4,
                "priority": "HIGH",
                "action": "INSTILL INITIAL ANTIDOTE THROUGH CATHETER (IF APPLICABLE)",
                "details": f"Instill portion of {primary_antidote.value.upper()} through the catheter lumen before removal.",
            })

        actions.extend([
            {
                "step_number": len(actions) + 1,
                "priority": "HIGH",
                "action": "REMOVE CATHETER & APPLY GENTLE PRESSURE",
                "details": "Gently withdraw the catheter. Apply light pressure with sterile gauze. Do NOT apply heavy friction or massage.",
            },
            {
                "step_number": len(actions) + 2,
                "priority": "HIGH",
                "action": "MARK EXTRAVASATION MARGINS",
                "details": "Outline the full visible area of erythema, induration, and edema with an indelible surgical skin marker for serial measurement.",
            },
            {
                "step_number": len(actions) + 3,
                "priority": "HIGH",
                "action": "INITIATE SPECIFIC ANTIDOTE PROTOCOL",
                "details": antidote_plan.dose_summary if antidote_plan else "No specific chemical antidote indicated; proceed with thermal management and symptomatic therapy.",
            },
            {
                "step_number": len(actions) + 4,
                "priority": "HIGH",
                "action": "APPLY THERMAL INTERVENTION",
                "details": f"{thermal_dict['protocol']}: {thermal_dict['frequency']}",
            },
            {
                "step_number": len(actions) + 5,
                "priority": "ROUTINE",
                "action": "ELEVATE AFFECTED EXTREMITY",
                "details": "Elevate the affected arm/limb above heart level for 48 hours to promote lymphatic drainage and minimize localized edema.",
            },
            {
                "step_number": len(actions) + 6,
                "priority": "URGENT" if ctcae["requires_urgent_surgical_review"] else "ROUTINE",
                "action": "SURGICAL CONSULTATION STATUS",
                "details": ctcae["surgical_consultation_status"],
            },
        ])

        # 5. Monitoring & Follow-up Schedule
        monitoring = [
            "Assess and document pain level (NRS 0-10), erythema diameter (cm), and edema every 15 minutes for the first 2 hours post-event.",
            "Inspect site every 4 hours for the subsequent 24 hours.",
            "Daily clinical examination and photographic documentation for Days 1 through 7.",
            "Outpatient clinical review at Day 14 and Day 28 post-extravasation.",
            "Instruct patient on red-flag signs (increasing pain, skin ulceration, numbness, paresthesias, fever > 38.0°C) requiring emergency presentation.",
        ]

        # 6. Documentation Requirements
        doc_reqs = [
            "Exact date, time, and clinical setting of extravasation detection.",
            f"Drug administered: {drug_name} (Estimated extravasated volume: {estimated_volume_ml} mL).",
            f"Vascular access device: {catheter_type} (Anatomical location and gauge).",
            "Patient symptoms reported (burning, stinging, pressure, numbness) and initial pain score.",
            "Exact interventions performed with timestamps (aspiration volume, antidote administered, compress applied).",
            "Baseline photographic recording with ruler / caliber scale next to marked perimeter.",
            "Institutional oncology incident report and pharmacy adverse drug event (ADE) reporting.",
        ]

        event_id = f"EXT-{datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%d%H%M%S')}"

        return ExtravasationEmergencyDossier(
            event_id=event_id,
            patient_id=patient_id,
            drug_name=drug_name,
            vesicant_class=vesicant_class,
            catheter_type=catheter_type,
            ctcae_severity=ctcae,
            thermal_protocol=thermal_dict,
            antidote_protocol=asdict(antidote_plan) if antidote_plan else None,
            ordered_action_checklist=actions,
            monitoring_and_followup_schedule=monitoring,
            documentation_requirements=doc_reqs,
        )


# ==============================================================================
# BATCH PROCESSING UTILITIES
# ==============================================================================

def process_batch_csv(input_csv_path: str, output_csv_path: str) -> int:
    """
    Reads a CSV of extravasation cases, assesses risk and management protocols,
    and writes out enriched clinical results.
    """
    engine = ChemotherapyExtravasationEngine()
    processed_count = 0

    with open(input_csv_path, mode="r", encoding="utf-8-sig") as infile:
        reader = csv.DictReader(infile)
        rows = list(reader)

    if not rows:
        return 0

    output_rows = []
    for row in rows:
        drug = row.get("drug", row.get("drug_name", "doxorubicin"))
        catheter = row.get("catheter", row.get("catheter_type", CatheterType.PERIPHERAL_FOREARM.value))
        vol = float(row.get("volume_ml", row.get("volume", 5.0)))
        elapsed = float(row.get("elapsed_hours", row.get("time_elapsed_hours", 0.5)))
        pain = int(row.get("pain_score", 4))
        ht = float(row.get("height_cm", 170.0))
        wt = float(row.get("weight_kg", 70.0))

        risk_eval = assess_extravasation_risk(drug, catheter)
        dossier = engine.evaluate_extravasation_event(
            drug_name=drug,
            catheter_type=catheter,
            estimated_volume_ml=vol,
            time_elapsed_hours=elapsed,
            patient_height_cm=ht,
            patient_weight_kg=wt,
            pain_score_0_to_10=pain,
        )

        out = dict(row)
        out["composite_risk_score"] = risk_eval.composite_risk_score
        out["risk_tier"] = risk_eval.risk_tier
        out["vesicant_class"] = dossier.vesicant_class
        out["ctcae_grade"] = dossier.ctcae_severity["ctcae_grade"]
        out["thermal_protocol"] = dossier.thermal_protocol["protocol"]
        out["antidote_indicated"] = dossier.antidote_protocol["antidote_name"] if dossier.antidote_protocol else "None"
        out["urgent_surgery_required"] = dossier.ctcae_severity["requires_urgent_surgical_review"]
        output_rows.append(out)
        processed_count += 1

    fieldnames = list(output_rows[0].keys())
    with open(output_csv_path, mode="w", encoding="utf-8", newline="") as outfile:
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(output_rows)

    return processed_count
