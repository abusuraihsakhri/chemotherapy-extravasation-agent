#!/usr/bin/env python3
"""
Unit Test Suite for Chemotherapy Extravasation Clinical Decision Support Engine.
"""

import os
import tempfile
import unittest
from chemotherapy_extravasation import (
    ChemotherapyExtravasationEngine,
    AntidoteCalculator,
    AntidoteType,
    CatheterType,
    VesicantClass,
    ThermalProtocol,
    DRUG_REGISTRY,
    assess_extravasation_risk,
    grade_ctcae_severity,
    calculate_bsa_mosteller,
    calculate_crcl_cockcroft_gault,
    process_batch_csv,
)


class TestDrugRegistry(unittest.TestCase):
    """Test drug classification and pharmacological registry."""

    def test_anthracycline_classification(self):
        dox = DRUG_REGISTRY["doxorubicin"]
        self.assertEqual(dox.vesicant_class, VesicantClass.DNA_BINDING_VESICANT)
        self.assertEqual(dox.primary_antidote, AntidoteType.DEXRAZOXANE)
        self.assertEqual(dox.thermal_protocol, ThermalProtocol.DRY_COLD)
        self.assertEqual(dox.antidote_window_hours, 6.0)

    def test_vinca_alkaloid_thermal_contraindication(self):
        vcr = DRUG_REGISTRY["vincristine"]
        self.assertEqual(vcr.vesicant_class, VesicantClass.NON_DNA_BINDING_VESICANT)
        self.assertEqual(vcr.primary_antidote, AntidoteType.HYALURONIDASE)
        self.assertEqual(vcr.thermal_protocol, ThermalProtocol.DRY_WARM)
        self.assertIn("COLD IS STRICTLY CONTRAINDICATED", vcr.thermal_rationale)

    def test_cisplatin_classification(self):
        cis = DRUG_REGISTRY["cisplatin"]
        self.assertEqual(cis.vesicant_class, VesicantClass.IRRITANT_WITH_VESICANT_POTENTIAL)
        self.assertEqual(cis.primary_antidote, AntidoteType.SODIUM_THIOSULFATE)
        self.assertEqual(cis.thermal_protocol, ThermalProtocol.DRY_COLD)

    def test_taxane_classification(self):
        pac = DRUG_REGISTRY["paclitaxel"]
        self.assertEqual(pac.vesicant_class, VesicantClass.NON_DNA_BINDING_VESICANT)
        self.assertEqual(pac.primary_antidote, AntidoteType.HYALURONIDASE)
        self.assertEqual(pac.thermal_protocol, ThermalProtocol.DRY_COLD)

    def test_oxaliplatin_thermal_rule(self):
        ox = DRUG_REGISTRY["oxaliplatin"]
        self.assertEqual(ox.thermal_protocol, ThermalProtocol.DRY_WARM)
        self.assertIn("cold exposure can trigger acute neurosensory dysesthesia", ox.thermal_rationale)


class TestPhysiologicalCalculations(unittest.TestCase):
    """Test BSA and Creatinine Clearance calculation accuracy."""

    def test_bsa_mosteller_standard(self):
        # 175 cm, 70 kg -> sqrt(175 * 70 / 3600) = sqrt(12250 / 3600) = sqrt(3.40277) = 1.8446 -> 1.84
        bsa = calculate_bsa_mosteller(175.0, 70.0)
        self.assertAlmostEqual(bsa, 1.84, places=2)

    def test_bsa_invalid_inputs(self):
        with self.assertRaises(ValueError):
            calculate_bsa_mosteller(-170, 70)
        with self.assertRaises(ValueError):
            calculate_bsa_mosteller(170, 0)

    def test_crcl_cockcroft_gault_male(self):
        # Male, Age 60, Wt 72 kg, Scr 1.0 -> ((140-60)*72) / (72*1.0) = 5760 / 72 = 80.0
        crcl = calculate_crcl_cockcroft_gault(60, 72.0, 1.0, is_female=False)
        self.assertEqual(crcl, 80.0)

    def test_crcl_cockcroft_gault_female(self):
        # Female, Age 60, Wt 72 kg, Scr 1.0 -> 80.0 * 0.85 = 68.0
        crcl = calculate_crcl_cockcroft_gault(60, 72.0, 1.0, is_female=True)
        self.assertEqual(crcl, 68.0)

    def test_crcl_invalid_inputs(self):
        with self.assertRaises(ValueError):
            calculate_crcl_cockcroft_gault(15, 70.0, 1.0, False)
        with self.assertRaises(ValueError):
            calculate_crcl_cockcroft_gault(60, -70.0, 1.0, False)


class TestAntidoteCalculators(unittest.TestCase):
    """Test antidote dosing logic, caps, and renal adjustments."""

    def setUp(self):
        self.calc = AntidoteCalculator()

    def test_dexrazoxane_standard_dosing(self):
        # BSA 1.8 m², CrCl 80 -> Day 1: 1800 mg, Day 2: 1800 mg, Day 3: 900 mg
        plan = self.calc.calculate_dexrazoxane(bsa_m2=1.8, crcl_ml_min=80.0, time_elapsed_hours=1.0)
        self.assertTrue(plan.is_within_window)
        self.assertFalse(plan.renal_adjustment_applied)
        self.assertEqual(plan.schedule[0]["dose_mg"], 1800)
        self.assertEqual(plan.schedule[1]["dose_mg"], 1800)
        self.assertEqual(plan.schedule[2]["dose_mg"], 900)

    def test_dexrazoxane_capping(self):
        # Large BSA 2.5 m² -> Day 1 uncapped = 2500 mg, capped at 2000 mg
        plan = self.calc.calculate_dexrazoxane(bsa_m2=2.5, crcl_ml_min=90.0, time_elapsed_hours=0.5)
        self.assertEqual(plan.schedule[0]["dose_mg"], 2000)
        self.assertEqual(plan.schedule[1]["dose_mg"], 2000)
        self.assertEqual(plan.schedule[2]["dose_mg"], 1000)

    def test_dexrazoxane_renal_adjustment(self):
        # CrCl 35 mL/min (< 50) -> 50% dose reduction
        # BSA 1.8 m² -> Day 1: 500 * 1.8 = 900 mg, Day 2: 900 mg, Day 3: 450 mg
        plan = self.calc.calculate_dexrazoxane(bsa_m2=1.8, crcl_ml_min=35.0, time_elapsed_hours=1.0)
        self.assertTrue(plan.renal_adjustment_applied)
        self.assertEqual(plan.schedule[0]["dose_mg"], 900)
        self.assertEqual(plan.schedule[1]["dose_mg"], 900)
        self.assertEqual(plan.schedule[2]["dose_mg"], 450)
        self.assertTrue(any("Renal impairment detected" in w for w in plan.contraindications_and_warnings))

    def test_dexrazoxane_expired_window_warning(self):
        plan = self.calc.calculate_dexrazoxane(bsa_m2=1.8, time_elapsed_hours=8.5)
        self.assertFalse(plan.is_within_window)
        self.assertTrue(any("exceeds the validated 6-hour" in w for w in plan.contraindications_and_warnings))

    def test_hyaluronidase_small_volume(self):
        plan = self.calc.calculate_hyaluronidase(estimated_volume_ml=1.5, time_elapsed_hours=0.5)
        self.assertEqual(plan.schedule[0]["total_units"], 150)
        self.assertEqual(plan.schedule[0]["num_sites"], 4)
        self.assertAlmostEqual(plan.schedule[0]["units_per_site"], 37.5)

    def test_hyaluronidase_moderate_volume(self):
        plan = self.calc.calculate_hyaluronidase(estimated_volume_ml=6.0, time_elapsed_hours=0.2)
        self.assertEqual(plan.schedule[0]["total_units"], 300)
        self.assertEqual(plan.schedule[0]["num_sites"], 6)
        self.assertEqual(plan.schedule[0]["units_per_site"], 50.0)

    def test_hyaluronidase_large_volume(self):
        plan = self.calc.calculate_hyaluronidase(estimated_volume_ml=25.0, time_elapsed_hours=0.5)
        self.assertEqual(plan.schedule[0]["total_units"], 750)
        self.assertEqual(plan.schedule[0]["num_sites"], 8)

    def test_sodium_thiosulfate_calculation(self):
        # 4.0 mL extravasate -> 2 mL/mL = 8.0 mL total STS
        plan = self.calc.calculate_sodium_thiosulfate(estimated_volume_ml=4.0, time_elapsed_hours=0.5)
        self.assertEqual(plan.schedule[0]["total_volume_ml"], 8.0)
        self.assertEqual(plan.schedule[0]["num_sites"], 5)
        self.assertEqual(plan.schedule[0]["volume_per_site_ml"], 1.6)

    def test_sodium_thiosulfate_capped_at_10ml(self):
        plan = self.calc.calculate_sodium_thiosulfate(estimated_volume_ml=12.0)
        self.assertEqual(plan.schedule[0]["total_volume_ml"], 10.0)

    def test_dmso_calculation(self):
        plan = self.calc.calculate_dmso(surface_area_cm2=30.0, time_elapsed_hours=1.0)
        self.assertEqual(plan.schedule[0]["dosage_drops"], 12)
        self.assertTrue(plan.is_within_window)


class TestCTCAESeverityStaging(unittest.TestCase):
    """Test CTCAE v5.0 grading and surgical consult determinations."""

    def test_grade_1_mild(self):
        res = grade_ctcae_severity(pain_score_0_to_10=2, erythema_present=True, edema_present=False)
        self.assertEqual(res["ctcae_grade"], 1)
        self.assertFalse(res["requires_urgent_surgical_review"])
        self.assertIn("ROUTINE", res["surgical_consultation_status"])

    def test_grade_2_moderate_blisters(self):
        res = grade_ctcae_severity(pain_score_0_to_10=5, blistering_present=True, blister_size_cm=0.5)
        self.assertEqual(res["ctcae_grade"], 2)
        self.assertFalse(res["requires_urgent_surgical_review"])
        self.assertIn("ADVISORY", res["surgical_consultation_status"])

    def test_grade_3_ulceration(self):
        res = grade_ctcae_severity(pain_score_0_to_10=8, ulceration_or_necrosis_present=True)
        self.assertEqual(res["ctcae_grade"], 3)
        self.assertTrue(res["requires_urgent_surgical_review"])
        self.assertIn("URGENT", res["surgical_consultation_status"])

    def test_grade_4_compartment_syndrome(self):
        res = grade_ctcae_severity(pain_score_0_to_10=9, compartment_syndrome_signs=True)
        self.assertEqual(res["ctcae_grade"], 4)
        self.assertTrue(res["requires_urgent_surgical_review"])
        self.assertIn("EMERGENT", res["surgical_consultation_status"])


class TestRiskAssessment(unittest.TestCase):
    """Test multifactorial extravasation risk scorer."""

    def test_low_risk_central_line_non_vesicant(self):
        eval_res = assess_extravasation_risk(
            drug_name="bleomycin",
            catheter_type=CatheterType.IMPLANTED_PORT.value,
        )
        self.assertEqual(eval_res.risk_tier, "LOW")
        self.assertLess(eval_res.composite_risk_score, 25.0)

    def test_high_risk_hand_vincristine(self):
        eval_res = assess_extravasation_risk(
            drug_name="vincristine",
            catheter_type=CatheterType.PERIPHERAL_HAND_WRIST.value,
            is_elderly_or_fragile_veins=True,
            multiple_venipuncture_attempts=True,
        )
        self.assertIn(eval_res.risk_tier, ("HIGH", "CRITICAL"))
        self.assertGreaterEqual(eval_res.composite_risk_score, 45.0)
        self.assertTrue(len(eval_res.clinical_flags) >= 2)

    def test_critical_risk_hand_doxorubicin_all_factors(self):
        eval_res = assess_extravasation_risk(
            drug_name="doxorubicin",
            catheter_type=CatheterType.PERIPHERAL_HAND_WRIST.value,
            multiple_venipuncture_attempts=True,
            is_elderly_or_fragile_veins=True,
            has_sensory_neuropathy=True,
            prior_radiation_to_limb=True,
            has_lymphedema=True,
            is_agitated_or_confused=True,
            is_infusion_pump_high_pressure=True,
            prolonged_infusion_gt_4h=True,
        )
        self.assertEqual(eval_res.risk_tier, "CRITICAL")
        self.assertGreaterEqual(eval_res.composite_risk_score, 70.0)


class TestFullEngineOrchestration(unittest.TestCase):
    """Test end-to-end extravasation clinical dossier creation."""

    def setUp(self):
        self.engine = ChemotherapyExtravasationEngine()

    def test_doxorubicin_extravasation_assessment(self):
        dossier = self.engine.evaluate_extravasation_event(
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
        self.assertEqual(dossier.vesicant_class, VesicantClass.DNA_BINDING_VESICANT.value)
        self.assertEqual(dossier.thermal_protocol["protocol"], "DRY COLD COMPRESS")
        self.assertIsNotNone(dossier.antidote_protocol)
        self.assertEqual(dossier.antidote_protocol["antidote_name"], "Dexrazoxane (Totect / Savene)")
        self.assertTrue(len(dossier.ordered_action_checklist) >= 6)
        self.assertEqual(dossier.ordered_action_checklist[0]["action"], "STOP INFUSION IMMEDIATELY")

    def test_vincristine_extravasation_warm_compress(self):
        dossier = self.engine.evaluate_extravasation_event(
            drug_name="vincristine",
            catheter_type="peripheral_hand_wrist",
            estimated_volume_ml=3.0,
            time_elapsed_hours=0.25,
            pain_score_0_to_10=4,
        )
        self.assertEqual(dossier.vesicant_class, VesicantClass.NON_DNA_BINDING_VESICANT.value)
        self.assertEqual(dossier.thermal_protocol["protocol"], "DRY WARM COMPRESS")
        self.assertIn("STRICTLY AVOID COLD", dossier.thermal_protocol["warning"])
        self.assertEqual(dossier.antidote_protocol["antidote_name"], "Hyaluronidase (Amphadase / Vitrase / Hylenex)")

    def test_batch_csv_processing(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            in_csv = os.path.join(tmpdir, "test_input.csv")
            out_csv = os.path.join(tmpdir, "test_output.csv")

            with open(in_csv, "w", encoding="utf-8") as f:
                f.write("drug,catheter,volume_ml,elapsed_hours,pain_score\n")
                f.write("doxorubicin,peripheral_hand_wrist,5.0,1.0,5\n")
                f.write("vincristine,peripheral_forearm,2.0,0.5,3\n")
                f.write("cisplatin,midline,15.0,1.5,4\n")

            count = process_batch_csv(in_csv, out_csv)
            self.assertEqual(count, 3)
            self.assertTrue(os.path.exists(out_csv))

            with open(out_csv, "r", encoding="utf-8") as f:
                lines = f.readlines()
                self.assertEqual(len(lines), 4)  # header + 3 rows


if __name__ == "__main__":
    unittest.main()
