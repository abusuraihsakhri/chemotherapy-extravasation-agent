#!/usr/bin/env python3
"""Regression tests for the chemotherapy extravasation reference tool."""

import csv
import os
import tempfile
import unittest

from chemotherapy_extravasation import (
    AntidoteCalculator,
    AntidoteType,
    CatheterType,
    ChemotherapyExtravasationEngine,
    DRUG_REGISTRY,
    ThermalProtocol,
    VesicantClass,
    assess_extravasation_risk,
    calculate_bsa_mosteller,
    calculate_crcl_cockcroft_gault,
    grade_ctcae_severity,
    process_batch_csv,
)


class TestDrugRegistry(unittest.TestCase):
    def test_anthracycline_classification(self):
        dox = DRUG_REGISTRY["doxorubicin"]
        self.assertEqual(dox.vesicant_class, VesicantClass.DNA_BINDING_VESICANT)
        self.assertEqual(dox.primary_antidote, AntidoteType.DEXRAZOXANE)
        self.assertEqual(dox.thermal_protocol, ThermalProtocol.DRY_COLD)

    def test_vinca_uses_hyaluronidase_and_warm_compress(self):
        vcr = DRUG_REGISTRY["vincristine"]
        self.assertEqual(vcr.vesicant_class, VesicantClass.NON_DNA_BINDING_VESICANT)
        self.assertEqual(vcr.primary_antidote, AntidoteType.HYALURONIDASE)
        self.assertEqual(vcr.thermal_protocol, ThermalProtocol.DRY_WARM)

    def test_taxane_is_irritant_with_vesicant_properties(self):
        pac = DRUG_REGISTRY["paclitaxel"]
        self.assertEqual(
            pac.vesicant_class, VesicantClass.IRRITANT_WITH_VESICANT_POTENTIAL
        )
        self.assertEqual(pac.primary_antidote, AntidoteType.HYALURONIDASE)
        self.assertEqual(pac.thermal_protocol, ThermalProtocol.DRY_WARM)

    def test_mitomycin_uses_dmso_pathway(self):
        mitomycin = DRUG_REGISTRY["mitomycin"]
        self.assertEqual(mitomycin.primary_antidote, AntidoteType.DMSO)


class TestPhysiologicalCalculations(unittest.TestCase):
    def test_bsa_mosteller(self):
        self.assertAlmostEqual(calculate_bsa_mosteller(175.0, 70.0), 1.84, places=2)

    def test_bsa_rejects_nonpositive_values(self):
        with self.assertRaises(ValueError):
            calculate_bsa_mosteller(-170, 70)
        with self.assertRaises(ValueError):
            calculate_bsa_mosteller(170, 0)

    def test_cockcroft_gault(self):
        self.assertEqual(
            calculate_crcl_cockcroft_gault(60, 72.0, 1.0, is_female=False), 80.0
        )
        self.assertEqual(
            calculate_crcl_cockcroft_gault(60, 72.0, 1.0, is_female=True), 68.0
        )

    def test_cockcroft_gault_rejects_invalid_inputs(self):
        with self.assertRaises(ValueError):
            calculate_crcl_cockcroft_gault(15, 70.0, 1.0, False)
        with self.assertRaises(ValueError):
            calculate_crcl_cockcroft_gault(60, -70.0, 1.0, False)


class TestAntidoteCalculators(unittest.TestCase):
    def setUp(self):
        self.calc = AntidoteCalculator()

    def test_dexrazoxane_standard_and_caps(self):
        plan = self.calc.calculate_dexrazoxane(
            bsa_m2=1.8, crcl_ml_min=80.0, time_elapsed_hours=1.0
        )
        self.assertEqual([d["dose_mg"] for d in plan.schedule], [1800, 1800, 900])
        capped = self.calc.calculate_dexrazoxane(
            bsa_m2=2.5, crcl_ml_min=90.0, time_elapsed_hours=0.5
        )
        self.assertEqual([d["dose_mg"] for d in capped.schedule], [2000, 2000, 1000])

    def test_dexrazoxane_renal_threshold_is_below_40(self):
        not_reduced = self.calc.calculate_dexrazoxane(
            bsa_m2=1.8, crcl_ml_min=40.0, time_elapsed_hours=1.0
        )
        reduced = self.calc.calculate_dexrazoxane(
            bsa_m2=1.8, crcl_ml_min=39.0, time_elapsed_hours=1.0
        )
        self.assertFalse(not_reduced.renal_adjustment_applied)
        self.assertTrue(reduced.renal_adjustment_applied)
        self.assertEqual(reduced.schedule[0]["dose_mg"], 900)

    def test_dexrazoxane_six_hour_window(self):
        plan = self.calc.calculate_dexrazoxane(bsa_m2=1.8, time_elapsed_hours=8.5)
        self.assertFalse(plan.is_within_window)
        self.assertTrue(
            any("6-hour" in w for w in plan.contraindications_and_warnings)
        )

    def test_hyaluronidase_is_not_volume_scaled(self):
        small = self.calc.calculate_hyaluronidase(
            estimated_volume_ml=1.5, time_elapsed_hours=0.5
        )
        large = self.calc.calculate_hyaluronidase(
            estimated_volume_ml=25.0, time_elapsed_hours=0.5
        )
        for plan in (small, large):
            self.assertEqual(plan.schedule[0]["total_units"], 150)
            self.assertEqual(plan.schedule[0]["num_sites"], 5)
            self.assertEqual(plan.schedule[0]["units_per_site"], 30.0)
            self.assertEqual(plan.schedule[0]["volume_per_site_ml"], 0.2)

    def test_sodium_thiosulfate_regimen(self):
        plan = self.calc.calculate_sodium_thiosulfate(
            estimated_volume_ml=25.0, time_elapsed_hours=0.5
        )
        self.assertEqual(plan.schedule[0]["total_volume_ml"], 1.0)
        self.assertEqual(plan.schedule[0]["num_sites"], 10)
        self.assertEqual(plan.schedule[0]["volume_per_site_ml"], 0.1)
        self.assertEqual(plan.schedule[0]["solution_concentration"], "1/6 M")

    def test_dmso_uses_measured_area(self):
        plan = self.calc.calculate_dmso(
            surface_area_cm2=30.0, time_elapsed_hours=0.2
        )
        self.assertEqual(plan.schedule[0]["area_coverage_cm2"], 60.0)
        self.assertTrue(plan.is_within_window)
        with self.assertRaises(ValueError):
            self.calc.calculate_dmso(surface_area_cm2=0)


class TestCTCAESeverityMapping(unittest.TestCase):
    def test_grade_1_painless_edema(self):
        res = grade_ctcae_severity(
            pain_score_0_to_10=0,
            erythema_present=False,
            edema_present=True,
        )
        self.assertEqual(res["ctcae_grade"], 1)
        self.assertFalse(res["requires_urgent_surgical_review"])

    def test_grade_2_erythema_with_symptoms(self):
        res = grade_ctcae_severity(
            pain_score_0_to_10=5,
            erythema_present=True,
            edema_present=True,
        )
        self.assertEqual(res["ctcae_grade"], 2)
        self.assertFalse(res["requires_urgent_surgical_review"])

    def test_grade_3_ulceration(self):
        res = grade_ctcae_severity(
            pain_score_0_to_10=3, ulceration_or_necrosis_present=True
        )
        self.assertEqual(res["ctcae_grade"], 3)
        self.assertTrue(res["requires_urgent_surgical_review"])

    def test_grade_4_life_threatening_findings(self):
        res = grade_ctcae_severity(
            pain_score_0_to_10=9, compartment_syndrome_signs=True
        )
        self.assertEqual(res["ctcae_grade"], 4)
        self.assertTrue(res["requires_urgent_surgical_review"])

    def test_invalid_pain_score_is_rejected(self):
        with self.assertRaises(ValueError):
            grade_ctcae_severity(pain_score_0_to_10=11)


class TestLegacyRiskHeuristic(unittest.TestCase):
    def test_known_drug_is_scored(self):
        result = assess_extravasation_risk(
            drug_name="bleomycin",
            catheter_type=CatheterType.IMPLANTED_PORT.value,
        )
        self.assertEqual(result.risk_tier, "LOW")

    def test_unknown_drug_is_not_inferred(self):
        with self.assertRaises(ValueError):
            assess_extravasation_risk(
                drug_name="not-a-real-agent",
                catheter_type=CatheterType.PERIPHERAL_FOREARM.value,
            )


class TestFullEngineOrchestration(unittest.TestCase):
    def setUp(self):
        self.engine = ChemotherapyExtravasationEngine()

    def test_unknown_drug_is_rejected(self):
        with self.assertRaises(ValueError):
            self.engine.evaluate_extravasation_event(
                drug_name="not-a-real-agent",
                catheter_type="peripheral_forearm",
                estimated_volume_ml=5.0,
                time_elapsed_hours=1.0,
            )

    def test_doxorubicin_with_measurements_calculates_dexrazoxane(self):
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
        self.assertEqual(
            dossier.vesicant_class, VesicantClass.DNA_BINDING_VESICANT.value
        )
        self.assertEqual(dossier.thermal_protocol["protocol"], "DRY COLD COMPRESS")
        self.assertTrue(dossier.antidote_protocol["is_indicated"])
        self.assertEqual(
            dossier.antidote_protocol["antidote_name"],
            "Dexrazoxane (Totect / Savene)",
        )
        self.assertGreater(len(dossier.antidote_protocol["schedule"]), 0)

    def test_doxorubicin_without_measurements_does_not_assume_bsa(self):
        dossier = self.engine.evaluate_extravasation_event(
            drug_name="doxorubicin",
            catheter_type="peripheral_forearm",
            estimated_volume_ml=5.0,
            time_elapsed_hours=1.0,
        )
        self.assertTrue(dossier.antidote_protocol["is_indicated"])
        self.assertEqual(dossier.antidote_protocol["schedule"], [])
        self.assertIn("height and weight", dossier.antidote_protocol["dose_summary"])

    def test_vincristine_uses_warm_compress_and_hyaluronidase(self):
        dossier = self.engine.evaluate_extravasation_event(
            drug_name="vincristine",
            catheter_type="peripheral_hand_wrist",
            estimated_volume_ml=3.0,
            time_elapsed_hours=0.25,
            pain_score_0_to_10=4,
        )
        self.assertEqual(dossier.thermal_protocol["protocol"], "DRY WARM COMPRESS")
        self.assertEqual(dossier.antidote_protocol["antidote_name"], "Hyaluronidase")

    def test_cisplatin_threshold_requires_both_volume_and_concentration(self):
        indicated = self.engine.evaluate_extravasation_event(
            drug_name="cisplatin",
            catheter_type="peripheral_forearm",
            estimated_volume_ml=25.0,
            time_elapsed_hours=0.5,
            drug_concentration_mg_ml=0.6,
        )
        self.assertEqual(
            indicated.vesicant_class, VesicantClass.DNA_BINDING_VESICANT.value
        )
        self.assertTrue(indicated.antidote_protocol["is_indicated"])
        self.assertEqual(
            indicated.antidote_protocol["antidote_name"],
            "Sodium Thiosulfate (1/6 M)",
        )

        low_volume = self.engine.evaluate_extravasation_event(
            drug_name="cisplatin",
            catheter_type="peripheral_forearm",
            estimated_volume_ml=20.0,
            time_elapsed_hours=0.5,
            drug_concentration_mg_ml=0.6,
        )
        self.assertFalse(low_volume.antidote_protocol["is_indicated"])

        low_concentration = self.engine.evaluate_extravasation_event(
            drug_name="cisplatin",
            catheter_type="peripheral_forearm",
            estimated_volume_ml=25.0,
            time_elapsed_hours=0.5,
            drug_concentration_mg_ml=0.4,
        )
        self.assertEqual(
            low_concentration.vesicant_class, VesicantClass.IRRITANT.value
        )
        self.assertFalse(low_concentration.antidote_protocol["is_indicated"])

    def test_batch_does_not_invent_patient_measurements(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            in_csv = os.path.join(tmpdir, "input.csv")
            out_csv = os.path.join(tmpdir, "output.csv")
            with open(in_csv, "w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(
                    ["drug", "catheter", "volume_ml", "elapsed_hours", "pain_score"]
                )
                writer.writerow(
                    ["doxorubicin", "peripheral_forearm", "5", "1", "5"]
                )

            count = process_batch_csv(in_csv, out_csv)
            self.assertEqual(count, 1)
            with open(out_csv, "r", encoding="utf-8", newline="") as handle:
                row = next(csv.DictReader(handle))
            self.assertEqual(row["antidote"], "Dexrazoxane")
            self.assertEqual(row["antidote_indicated"], "True")


if __name__ == "__main__":
    unittest.main()
