#!/usr/bin/env python3
"""
Command-Line Interface for Chemotherapy Extravasation Clinical Decision Support Engine.

Usage:
    python cli.py assess --drug doxorubicin --catheter peripheral_forearm --volume 8.0 --elapsed 1.5 --weight 72 --height 175 --age 58 --creatinine 1.1
    python cli.py antidote --type dexrazoxane --bsa 1.85 --crcl 45.0 --elapsed 2.0
    python cli.py antidote --type hyaluronidase --volume 6.0
    python cli.py risk --drug vincristine --catheter peripheral_hand_wrist --fragile-veins --multiple-attempts
    python cli.py stage --pain 8 --ulceration
    python cli.py batch --input sample.csv --output results.csv
    python cli.py interactive
"""

import argparse
import json
import os
import sys
from dataclasses import asdict

from chemotherapy_extravasation import (
    ChemotherapyExtravasationEngine,
    AntidoteCalculator,
    AntidoteType,
    CatheterType,
    VesicantClass,
    DRUG_REGISTRY,
    assess_extravasation_risk,
    grade_ctcae_severity,
    calculate_bsa_mosteller,
    calculate_crcl_cockcroft_gault,
    process_batch_csv,
)


def format_dossier_display(dossier):
    """Formats and prints an emergency clinical dossier."""
    d = asdict(dossier)
    print("=" * 80)
    print(f"  CHEMOTHERAPY EXTRAVASATION EMERGENCY CLINICAL DOSSIER [{d['event_id']}]")
    print("=" * 80)
    print(f"  Patient ID:       {d['patient_id']}")
    print(f"  Drug Name:        {d['drug_name'].upper()}")
    print(f"  Vesicant Class:   {d['vesicant_class'].upper().replace('_', ' ')}")
    print(f"  Catheter Type:    {d['catheter_type']}")
    print(f"  Timestamp (UTC):  {d['timestamp_utc']}")
    print("-" * 80)
    
    ctcae = d['ctcae_severity']
    print(f"  CTCAE v5.0 Grade: Grade {ctcae['ctcae_grade']}")
    print(f"  Severity Detail:  {ctcae['description']}")
    print(f"  Surgical Consult: {ctcae['surgical_consultation_status']}")
    print("-" * 80)

    thermal = d['thermal_protocol']
    print(f"  THERMAL INTERVENTION: {thermal['protocol']}")
    print(f"  Frequency/Regimen:    {thermal['frequency']}")
    print(f"  Clinical Rationale:   {thermal['rationale']}")
    print(f"  Warning / Note:       {thermal['warning']}")
    print("-" * 80)

    antidote = d.get('antidote_protocol')
    if antidote:
        print(f"  SPECIFIC ANTIDOTE:    {antidote['antidote_name']}")
        print(f"  Indication Status:    {'INDICATED' if antidote['is_indicated'] else 'NOT INDICATED'}")
        if antidote['urgency_window_hours'] > 0:
            print(f"  Therapeutic Window:   Within {antidote['urgency_window_hours']} hours (Elapsed: {antidote['time_elapsed_hours']}h - {'OK' if antidote['is_within_window'] else 'OUTSIDE'})")
        else:
            print("  Therapeutic Window:   Verify timing in the agent-specific institutional protocol")
        print(f"  Dosing Summary:       {antidote['dose_summary']}")
        if antidote.get('renal_adjustment_applied'):
            print("  Renal Adjustment:     APPLIED (Dose reduced by 50% for CrCl < 40 mL/min)")
        print("\n  Administration Steps:")
        for instr in antidote['administration_instructions']:
            print(f"    - {instr}")
        if antidote.get('contraindications_and_warnings'):
            print("\n  Antidote Warnings:")
            for w in antidote['contraindications_and_warnings']:
                print(f"    ! {w}")
    else:
        print("  SPECIFIC ANTIDOTE:    None Indicated (Supportive care and thermal intervention only)")
    print("-" * 80)

    print("  PRIORITIZED ACTION CHECKLIST:")
    for step in d['ordered_action_checklist']:
        print(f"  [{step['step_number']}] [{step['priority']}] {step['action']}")
        print(f"      {step['details']}")
    print("-" * 80)

    print("  MONITORING & FOLLOW-UP SCHEDULE:")
    for m in d['monitoring_and_followup_schedule']:
        print(f"    * {m}")
    print("=" * 80)


def cmd_assess(args):
    engine = ChemotherapyExtravasationEngine()
    dossier = engine.evaluate_extravasation_event(
        drug_name=args.drug,
        catheter_type=args.catheter,
        estimated_volume_ml=args.volume,
        time_elapsed_hours=args.elapsed,
        drug_concentration_mg_ml=args.concentration,
        extravasation_surface_area_cm2=args.surface_area,
        patient_height_cm=args.height,
        patient_weight_kg=args.weight,
        patient_age_years=args.age,
        serum_creatinine_mg_dl=args.creatinine,
        is_female=args.female,
        pain_score_0_to_10=args.pain,
        erythema_present=args.erythema,
        edema_present=args.edema,
        blistering_present=args.blistering,
        ulceration_or_necrosis_present=args.ulceration,
        tissue_sloughing_or_eschar=args.sloughing,
        compartment_syndrome_signs=args.compartment,
        loss_of_extremity_function=args.loss_of_function,
        patient_id=args.patient_id or "PATIENT-001",
    )

    if args.json:
        print(json.dumps(asdict(dossier), indent=2, default=str))
    else:
        format_dossier_display(dossier)
    return 0


def cmd_antidote(args):
    calc = AntidoteCalculator()
    antidote_type = args.type.lower()
    
    if antidote_type in ("dexrazoxane", "totect", "savene"):
        bsa = args.bsa
        if bsa is None and args.height is not None and args.weight is not None:
            bsa = calculate_bsa_mosteller(args.height, args.weight)
        elif bsa is None:
            raise ValueError("Dexrazoxane dosing requires --bsa or both --height and --weight.")

        crcl = args.crcl
        if crcl is None and args.age and args.weight and args.creatinine:
            crcl = calculate_crcl_cockcroft_gault(args.age, args.weight, args.creatinine, args.female)

        plan = calc.calculate_dexrazoxane(bsa_m2=bsa, crcl_ml_min=crcl, time_elapsed_hours=args.elapsed or 0.0)

    elif antidote_type in ("hyaluronidase", "vitrase", "amphadase", "hylenex"):
        plan = calc.calculate_hyaluronidase(estimated_volume_ml=args.volume or 5.0, time_elapsed_hours=args.elapsed or 0.0)

    elif antidote_type in ("sodium_thiosulfate", "sts", "thiosulfate"):
        plan = calc.calculate_sodium_thiosulfate(estimated_volume_ml=args.volume or 5.0, time_elapsed_hours=args.elapsed or 0.0)

    elif antidote_type in ("dmso", "dimethyl_sulfoxide"):
        if args.surface_area is None:
            raise ValueError("DMSO calculation requires --surface-area in cm².")
        plan = calc.calculate_dmso(surface_area_cm2=args.surface_area, time_elapsed_hours=args.elapsed or 0.0)

    else:
        print(f"Error: Unknown antidote type '{args.type}'. Supported: dexrazoxane, hyaluronidase, sodium_thiosulfate, dmso", file=sys.stderr)
        return 1

    if args.json:
        print(json.dumps(asdict(plan), indent=2, default=str))
    else:
        print("=" * 70)
        print(f"  ANTIDOTE PROTOCOL: {plan.antidote_name.upper()}")
        print("=" * 70)
        print(f"  Dosing Summary:     {plan.dose_summary}")
        print(f"  Therapeutic Window: {plan.urgency_window_hours} hours (Elapsed: {plan.time_elapsed_hours}h - {'VALID' if plan.is_within_window else 'EXPIRED'})")
        print(f"  Renal Adjustment:   {'YES (50% reduction)' if plan.renal_adjustment_applied else 'No'}")
        print("\n  Instructions:")
        for instr in plan.administration_instructions:
            print(f"    - {instr}")
        if plan.contraindications_and_warnings:
            print("\n  Warnings & Precautions:")
            for w in plan.contraindications_and_warnings:
                print(f"    ! {w}")
        print("=" * 70)
    return 0


def cmd_risk(args):
    result = assess_extravasation_risk(
        drug_name=args.drug,
        catheter_type=args.catheter,
        multiple_venipuncture_attempts=args.multiple_attempts,
        is_elderly_or_fragile_veins=args.fragile_veins,
        has_sensory_neuropathy=args.neuropathy,
        prior_radiation_to_limb=args.prior_radiation,
        has_lymphedema=args.lymphedema,
        is_agitated_or_confused=args.agitated,
        is_infusion_pump_high_pressure=args.pump_pressure,
        prolonged_infusion_gt_4h=args.prolonged,
    )

    if args.json:
        print(json.dumps(asdict(result), indent=2, default=str))
    else:
        print("=" * 70)
        print("  UNVALIDATED EXTRAVASATION RISK HEURISTIC")
        print("=" * 70)
        print(f"  Drug Name:            {result.drug_name.upper()} ({result.vesicant_class})")
        print(f"  Catheter Type:        {result.catheter_type}")
        print(f"  Composite Risk Score: {result.composite_risk_score} / 100")
        print(f"  Heuristic Tier:       [{result.risk_tier}]")
        print("  NOTE: This is not a validated clinical prediction rule and must not determine treatment or vascular access.")
        print("\n  Heuristic Breakdown:")
        for k, v in result.score_breakdown.items():
            print(f"    - {k.replace('_', ' ').title()}: {v} pts")
        if result.clinical_flags:
            print("\n  Clinical Risk Flags:")
            for flag in result.clinical_flags:
                print(f"    ! {flag}")
        print("\n  Prevention Recommendations:")
        for rec in result.prevention_recommendations:
            print(f"    * {rec}")
        print("=" * 70)
    return 0


def cmd_stage(args):
    result = grade_ctcae_severity(
        pain_score_0_to_10=args.pain,
        erythema_present=args.erythema,
        edema_present=args.edema,
        blistering_present=args.blistering,
        blister_size_cm=args.blister_size,
        ulceration_or_necrosis_present=args.ulceration,
        tissue_sloughing_or_eschar=args.sloughing,
        compartment_syndrome_signs=args.compartment,
        loss_of_extremity_function=args.loss_of_function,
    )

    if args.json:
        print(json.dumps(result, indent=2, default=str))
    else:
        print("=" * 70)
        print("  CTCAE v5.0 EXTRAVASATION SEVERITY STAGING")
        print("=" * 70)
        print(f"  Grade:                    Grade {result['ctcae_grade']}")
        print(f"  Description:              {result['description']}")
        print(f"  Surgical Action:          {result['surgical_consultation_status']}")
        print(f"  Urgent Surgical Review:   {'YES' if result['requires_urgent_surgical_review'] else 'No'}")
        print("=" * 70)
    return 0


def cmd_batch(args):
    if not os.path.exists(args.input):
        print(f"Error: Input file '{args.input}' not found.", file=sys.stderr)
        return 1
    count = process_batch_csv(args.input, args.output)
    print(f"Successfully processed {count} extravasation records from '{args.input}' -> '{args.output}'.")
    return 0


def cmd_interactive(args):
    print("=" * 70)
    print("  CHEMOTHERAPY EXTRAVASATION INTERACTIVE DECISION SUPPORT")
    print("=" * 70)
    drug = input("Enter Chemotherapy Drug Name (e.g., doxorubicin, vincristine, cisplatin): ").strip()
    if not drug:
        print("Drug name is required.")
        return 1
    
    print("\nSelect Catheter Type:")
    catheters = [
        ("1", CatheterType.PERIPHERAL_HAND_WRIST.value),
        ("2", CatheterType.PERIPHERAL_FOREARM.value),
        ("3", CatheterType.PERIPHERAL_ANTECUBITAL.value),
        ("4", CatheterType.MIDLINE.value),
        ("5", CatheterType.PICC.value),
        ("6", CatheterType.IMPLANTED_PORT.value),
    ]
    for num, cat in catheters:
        print(f"  [{num}] {cat}")
    c_choice = input("Choice (default 2): ").strip() or "2"
    catheter_map = dict(catheters)
    catheter_val = catheter_map.get(c_choice, CatheterType.PERIPHERAL_FOREARM.value)

    vol_str = input("\nEstimated volume extravasated in mL (default 5.0): ").strip() or "5.0"
    vol = float(vol_str)

    elapsed_str = input("Time elapsed since extravasation in hours (default 0.5): ").strip() or "0.5"
    elapsed = float(elapsed_str)

    pain_str = input("Patient pain score (0-10, default 5): ").strip() or "5"
    pain = int(pain_str)

    wt_str = input("Patient weight in kg (optional): ").strip()
    wt = float(wt_str) if wt_str else None

    ht_str = input("Patient height in cm (optional): ").strip()
    ht = float(ht_str) if ht_str else None

    concentration_str = input("Drug concentration in mg/mL if relevant (optional): ").strip()
    concentration = float(concentration_str) if concentration_str else None

    area_str = input("Measured extravasation area in cm² if relevant (optional): ").strip()
    surface_area = float(area_str) if area_str else None

    blister = input("Are blisters present? (y/n, default n): ").strip().lower().startswith("y")
    ulcer = input("Is skin ulceration/necrosis present? (y/n, default n): ").strip().lower().startswith("y")

    engine = ChemotherapyExtravasationEngine()
    dossier = engine.evaluate_extravasation_event(
        drug_name=drug,
        catheter_type=catheter_val,
        estimated_volume_ml=vol,
        time_elapsed_hours=elapsed,
        drug_concentration_mg_ml=concentration,
        extravasation_surface_area_cm2=surface_area,
        patient_height_cm=ht,
        patient_weight_kg=wt,
        pain_score_0_to_10=pain,
        blistering_present=blister,
        ulceration_or_necrosis_present=ulcer,
    )
    print("\n")
    format_dossier_display(dossier)
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="chemotherapy-extravasation",
        description="Reference tool for antineoplastic extravasation management using current ONS/ASCO guidance and CTCAE v5.0 mapping.",
    )
    parser.add_argument("--json", action="store_true", help="Output results in JSON format")

    subparsers = parser.add_subparsers(dest="command")

    # Assess
    p_assess = subparsers.add_parser("assess", help="Comprehensive extravasation event assessment & protocol generation")
    p_assess.add_argument("--drug", required=True, help="Antineoplastic agent generic name")
    p_assess.add_argument("--catheter", default="peripheral_forearm", help="Catheter / access device type")
    p_assess.add_argument("--volume", type=float, default=5.0, help="Estimated volume extravasated in mL")
    p_assess.add_argument("--elapsed", type=float, default=0.5, help="Time elapsed since event in hours")
    p_assess.add_argument("--concentration", type=float, help="Drug concentration in mg/mL (required to determine cisplatin thresholds)")
    p_assess.add_argument("--surface-area", type=float, help="Measured extravasation surface area in cm² (required for DMSO calculation)")
    p_assess.add_argument("--pain", type=int, default=4, help="Patient pain score (0-10)")
    p_assess.add_argument("--erythema", action=argparse.BooleanOptionalAction, default=True, help="Erythema present (use --no-erythema when absent)")
    p_assess.add_argument("--edema", action=argparse.BooleanOptionalAction, default=True, help="Edema present (use --no-edema when absent)")
    p_assess.add_argument("--height", type=float, help="Patient height in cm")
    p_assess.add_argument("--weight", type=float, help="Patient weight in kg")
    p_assess.add_argument("--age", type=int, help="Patient age in years")
    p_assess.add_argument("--creatinine", type=float, help="Serum creatinine in mg/dL")
    p_assess.add_argument("--female", action="store_true", help="Apply the female coefficient in Cockcroft-Gault")
    p_assess.add_argument("--blistering", action="store_true", help="Blisters present at site")
    p_assess.add_argument("--ulceration", action="store_true", help="Ulceration / necrosis present")
    p_assess.add_argument("--sloughing", action="store_true", help="Tissue sloughing or eschar present")
    p_assess.add_argument("--compartment", action="store_true", help="Compartment syndrome signs present")
    p_assess.add_argument("--loss-of-function", action="store_true", help="Loss of extremity function present")
    p_assess.add_argument("--patient-id", help="Clinical patient identifier")

    # Antidote
    p_antidote = subparsers.add_parser("antidote", help="Calculate specific antidote dosing regimen")
    p_antidote.add_argument("--type", required=True, choices=["dexrazoxane", "hyaluronidase", "sodium_thiosulfate", "dmso"], help="Antidote type")
    p_antidote.add_argument("--bsa", type=float, help="Body surface area in m² (for dexrazoxane)")
    p_antidote.add_argument("--crcl", type=float, help="Creatinine clearance in mL/min")
    p_antidote.add_argument("--volume", type=float, default=5.0, help="Extravasated volume in mL")
    p_antidote.add_argument("--elapsed", type=float, default=0.5, help="Hours elapsed since event")
    p_antidote.add_argument("--height", type=float, help="Height in cm")
    p_antidote.add_argument("--weight", type=float, help="Weight in kg")
    p_antidote.add_argument("--age", type=int, help="Age in years")
    p_antidote.add_argument("--creatinine", type=float, help="Serum creatinine in mg/dL")
    p_antidote.add_argument("--female", action="store_true", help="Apply the female coefficient in Cockcroft-Gault")
    p_antidote.add_argument("--surface-area", type=float, help="Measured extravasation surface area in cm² (required for DMSO)")

    # Risk
    p_risk = subparsers.add_parser("risk", help="Legacy unvalidated pre-infusion risk heuristic")
    p_risk.add_argument("--drug", required=True, help="Chemotherapy drug name")
    p_risk.add_argument("--catheter", default="peripheral_forearm", help="Catheter type")
    p_risk.add_argument("--multiple-attempts", action="store_true", help="Multiple venipuncture attempts made")
    p_risk.add_argument("--fragile-veins", action="store_true", help="Fragile or sclerotic veins")
    p_risk.add_argument("--neuropathy", action="store_true", help="Sensory neuropathy present")
    p_risk.add_argument("--prior-radiation", action="store_true", help="Prior radiation to limb")
    p_risk.add_argument("--lymphedema", action="store_true", help="Lymphedema present")
    p_risk.add_argument("--agitated", action="store_true", help="Agitated or confused patient")
    p_risk.add_argument("--pump-pressure", action="store_true", help="High-pressure infusion pump")
    p_risk.add_argument("--prolonged", action="store_true", help="Prolonged infusion > 4 hours")

    # Stage
    p_stage = subparsers.add_parser("stage", help="Grade CTCAE v5.0 severity")
    p_stage.add_argument("--pain", type=int, default=3, help="Pain score (0-10)")
    p_stage.add_argument("--erythema", action=argparse.BooleanOptionalAction, default=True, help="Erythema present (use --no-erythema when absent)")
    p_stage.add_argument("--edema", action=argparse.BooleanOptionalAction, default=True, help="Edema present (use --no-edema when absent)")
    p_stage.add_argument("--blistering", action="store_true", help="Blistering present")
    p_stage.add_argument("--blister-size", type=float, default=0.0, help="Blister diameter in cm")
    p_stage.add_argument("--ulceration", action="store_true", help="Ulceration or necrosis present")
    p_stage.add_argument("--sloughing", action="store_true", help="Eschar / tissue sloughing")
    p_stage.add_argument("--compartment", action="store_true", help="Compartment syndrome signs")
    p_stage.add_argument("--loss-of-function", action="store_true", help="Loss of limb function")

    # Batch
    p_batch = subparsers.add_parser("batch", help="Process batch cases from CSV")
    p_batch.add_argument("-i", "--input", required=True, help="Input CSV file path")
    p_batch.add_argument("-o", "--output", default="extravasation_results.csv", help="Output CSV file path")

    # Interactive
    p_interactive = subparsers.add_parser("interactive", help="Interactive clinical wizard")

    args = parser.parse_args(argv)

    try:
        if args.command == "assess":
            return cmd_assess(args)
        if args.command == "antidote":
            return cmd_antidote(args)
        if args.command == "risk":
            return cmd_risk(args)
        if args.command == "stage":
            return cmd_stage(args)
        if args.command == "batch":
            return cmd_batch(args)
        if args.command == "interactive":
            return cmd_interactive(args)
        parser.print_help()
        return 0
    except (ValueError, OSError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
