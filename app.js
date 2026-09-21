"use strict";

const PYODIDE_INDEX = "https://cdn.jsdelivr.net/pyodide/v314.0.7/full/";

const elements = {
  form: document.getElementById("assessment-form"),
  drug: document.getElementById("drug"),
  catheter: document.getElementById("catheter"),
  volume: document.getElementById("volume"),
  elapsed: document.getElementById("elapsed"),
  pain: document.getElementById("pain"),
  erythema: document.getElementById("erythema"),
  edema: document.getElementById("edema"),
  blistering: document.getElementById("blistering"),
  ulceration: document.getElementById("ulceration"),
  sloughing: document.getElementById("sloughing"),
  compartment: document.getElementById("compartment"),
  lossFunction: document.getElementById("loss-function"),
  concentration: document.getElementById("concentration"),
  surfaceArea: document.getElementById("surface-area"),
  height: document.getElementById("height"),
  weight: document.getElementById("weight"),
  age: document.getElementById("age"),
  creatinine: document.getElementById("creatinine"),
  femaleCoefficient: document.getElementById("female-coefficient"),
  assessButton: document.getElementById("assess-button"),
  resetButton: document.getElementById("reset-button"),
  themeToggle: document.getElementById("theme-toggle"),
  runtimeStatus: document.getElementById("runtime-status"),
  formError: document.getElementById("form-error"),
  emptyState: document.getElementById("empty-state"),
  resultContent: document.getElementById("result-content"),
  eventId: document.getElementById("event-id"),
  urgentBanner: document.getElementById("urgent-banner"),
  ctcaeGrade: document.getElementById("ctcae-grade"),
  ctcaeDetail: document.getElementById("ctcae-detail"),
  classification: document.getElementById("classification"),
  thermalName: document.getElementById("thermal-name"),
  thermalFrequency: document.getElementById("thermal-frequency"),
  antidoteName: document.getElementById("antidote-name"),
  antidoteSummary: document.getElementById("antidote-summary"),
  thermalRationale: document.getElementById("thermal-rationale"),
  thermalWarning: document.getElementById("thermal-warning"),
  actionList: document.getElementById("action-list"),
  antidoteInstructions: document.getElementById("antidote-instructions"),
  antidoteWarnings: document.getElementById("antidote-warnings"),
  monitoringList: document.getElementById("monitoring-list"),
  documentationList: document.getElementById("documentation-list"),
};

let pyodide = null;

function setRuntimeStatus(message, state = "") {
  elements.runtimeStatus.textContent = message;
  elements.runtimeStatus.classList.remove("ready", "error");
  if (state) {
    elements.runtimeStatus.classList.add(state);
  }
}

function optionalNumber(input) {
  const value = input.value.trim();
  return value === "" ? null : Number(value);
}

function titleCase(value) {
  return String(value)
    .replaceAll("_", " ")
    .replace(/\b\w/g, (letter) => letter.toUpperCase());
}

function setTheme(theme) {
  document.documentElement.dataset.theme = theme;
  elements.themeToggle.textContent = theme === "dark" ? "Light" : "Dark";
  elements.themeToggle.setAttribute(
    "aria-label",
    theme === "dark" ? "Switch to light theme" : "Switch to dark theme"
  );
  localStorage.setItem("extravasation-theme", theme);
}

function initializeTheme() {
  const saved = localStorage.getItem("extravasation-theme");
  const preferred =
    window.matchMedia &&
    window.matchMedia("(prefers-color-scheme: dark)").matches
      ? "dark"
      : "light";
  setTheme(saved === "dark" || saved === "light" ? saved : preferred);
}

function appendListItems(container, items, fallback = "None generated.") {
  container.replaceChildren();
  const values = Array.isArray(items) && items.length ? items : [fallback];
  for (const value of values) {
    const li = document.createElement("li");
    li.textContent = String(value);
    container.append(li);
  }
}

function renderActions(actions) {
  elements.actionList.replaceChildren();
  for (const step of actions) {
    const li = document.createElement("li");

    const title = document.createElement("span");
    title.className = "action-title";
    title.textContent = step.action;

    const priority = document.createElement("span");
    priority.className = "action-priority";
    priority.textContent = step.priority;
    title.append(priority);

    const detail = document.createElement("span");
    detail.className = "action-detail";
    detail.textContent = step.details;

    li.append(title, detail);
    elements.actionList.append(li);
  }
}

function clearError() {
  elements.formError.hidden = true;
  elements.formError.textContent = "";
}

function showError(message) {
  elements.formError.textContent = message;
  elements.formError.hidden = false;
}

function clearResults() {
  elements.resultContent.hidden = true;
  elements.emptyState.hidden = false;
  elements.eventId.textContent = "";
  elements.urgentBanner.hidden = true;
  elements.urgentBanner.textContent = "";
}

function validatePayload(payload) {
  if (!payload.drug_name) {
    throw new Error("Select an antineoplastic agent.");
  }
  if (!(payload.estimated_volume_ml > 0)) {
    throw new Error("Estimated volume must be greater than 0 mL.");
  }
  if (!(payload.time_elapsed_hours >= 0)) {
    throw new Error("Elapsed time cannot be negative.");
  }
  if (
    !Number.isInteger(payload.pain_score_0_to_10) ||
    payload.pain_score_0_to_10 < 0 ||
    payload.pain_score_0_to_10 > 10
  ) {
    throw new Error("Pain score must be an integer from 0 to 10.");
  }
}

function buildPayload() {
  const payload = {
    drug_name: elements.drug.value,
    catheter_type: elements.catheter.value,
    estimated_volume_ml: Number(elements.volume.value),
    time_elapsed_hours: Number(elements.elapsed.value),
    drug_concentration_mg_ml: optionalNumber(elements.concentration),
    extravasation_surface_area_cm2: optionalNumber(elements.surfaceArea),
    patient_height_cm: optionalNumber(elements.height),
    patient_weight_kg: optionalNumber(elements.weight),
    patient_age_years: optionalNumber(elements.age),
    serum_creatinine_mg_dl: optionalNumber(elements.creatinine),
    is_female: elements.femaleCoefficient.checked,
    pain_score_0_to_10: Number(elements.pain.value),
    erythema_present: elements.erythema.checked,
    edema_present: elements.edema.checked,
    blistering_present: elements.blistering.checked,
    ulceration_or_necrosis_present: elements.ulceration.checked,
    tissue_sloughing_or_eschar: elements.sloughing.checked,
    compartment_syndrome_signs: elements.compartment.checked,
    loss_of_extremity_function: elements.lossFunction.checked,
    patient_id: "BROWSER-SESSION",
  };

  if (payload.patient_age_years !== null) {
    payload.patient_age_years = Number.parseInt(
      String(payload.patient_age_years),
      10
    );
  }
  payload.pain_score_0_to_10 = Number.parseInt(
    String(payload.pain_score_0_to_10),
    10
  );

  validatePayload(payload);
  return payload;
}

function renderResult(data) {
  elements.emptyState.hidden = true;
  elements.resultContent.hidden = false;
  elements.eventId.textContent = data.event_id || "";

  const ctcae = data.ctcae_severity;
  elements.ctcaeGrade.textContent = `Grade ${ctcae.ctcae_grade}`;
  elements.ctcaeDetail.textContent = ctcae.description;
  elements.classification.textContent = titleCase(data.vesicant_class);

  elements.thermalName.textContent = data.thermal_protocol.protocol;
  elements.thermalFrequency.textContent = data.thermal_protocol.frequency;
  elements.thermalRationale.textContent = data.thermal_protocol.rationale;
  elements.thermalWarning.textContent = data.thermal_protocol.warning;

  const antidote = data.antidote_protocol;
  if (antidote) {
    elements.antidoteName.textContent = antidote.is_indicated
      ? antidote.antidote_name
      : "Not indicated";
    elements.antidoteSummary.textContent = antidote.is_indicated
      ? antidote.dose_summary
      : `${antidote.antidote_name}: ${antidote.dose_summary}`;
    appendListItems(
      elements.antidoteInstructions,
      antidote.administration_instructions,
      "No administration instructions generated."
    );
    appendListItems(
      elements.antidoteWarnings,
      antidote.contraindications_and_warnings,
      "No additional antidote warnings generated."
    );
  } else {
    elements.antidoteName.textContent = "No specific antidote";
    elements.antidoteSummary.textContent =
      "Follow the agent-specific thermal and supportive-care pathway.";
    appendListItems(
      elements.antidoteInstructions,
      [],
      "No specific antidote administration is generated for this agent."
    );
    appendListItems(elements.antidoteWarnings, [], "No antidote warning.");
  }

  if (ctcae.requires_urgent_surgical_review) {
    elements.urgentBanner.textContent =
      ctcae.surgical_consultation_status;
    elements.urgentBanner.hidden = false;
  } else {
    elements.urgentBanner.hidden = true;
    elements.urgentBanner.textContent = "";
  }

  renderActions(data.ordered_action_checklist || []);
  appendListItems(elements.monitoringList, data.monitoring_and_followup_schedule);
  appendListItems(elements.documentationList, data.documentation_requirements);
}

async function runAssessment(payload) {
  pyodide.globals.set("payload_json", JSON.stringify(payload));
  const result = await pyodide.runPythonAsync(`
import json
from dataclasses import asdict
from chemotherapy_extravasation import ChemotherapyExtravasationEngine

_payload = json.loads(payload_json)
_result = ChemotherapyExtravasationEngine().evaluate_extravasation_event(**_payload)
json.dumps(asdict(_result))
`);
  return JSON.parse(result);
}

function populateDrugSelect(drugs) {
  elements.drug.replaceChildren();
  for (const drug of drugs) {
    const option = document.createElement("option");
    option.value = drug;
    option.textContent = titleCase(drug);
    if (drug === "doxorubicin") {
      option.defaultSelected = true;
      option.selected = true;
    }
    elements.drug.append(option);
  }
  elements.drug.disabled = false;
}

function shortError(error) {
  const raw = String(error && error.message ? error.message : error);
  const lines = raw
    .split("\n")
    .map((line) => line.trim())
    .filter(Boolean);
  return lines.at(-1) || "Assessment failed.";
}

async function initializePython() {
  try {
    setRuntimeStatus("Loading Python…");
    pyodide = await loadPyodide({ indexURL: PYODIDE_INDEX });

    const response = await fetch("./chemotherapy_extravasation.py", {
      cache: "no-store",
    });
    if (!response.ok) {
      throw new Error(
        `Unable to load clinical engine (HTTP ${response.status}).`
      );
    }

    const source = await response.text();
    pyodide.FS.writeFile("chemotherapy_extravasation.py", source);
    await pyodide.runPythonAsync("import chemotherapy_extravasation");

    const drugsJson = pyodide.runPython(`
import json
from chemotherapy_extravasation import DRUG_REGISTRY
json.dumps(sorted(DRUG_REGISTRY.keys()))
`);
    populateDrugSelect(JSON.parse(drugsJson));

    elements.assessButton.disabled = false;
    setRuntimeStatus("Python ready", "ready");
  } catch (error) {
    setRuntimeStatus("Runtime failed", "error");
    showError(
      `${shortError(error)} Reload the page after checking network access to the Pyodide CDN.`
    );
  }
}

elements.form.addEventListener("submit", async (event) => {
  event.preventDefault();
  clearError();

  if (!pyodide) {
    showError("Python runtime is not ready.");
    return;
  }

  try {
    const payload = buildPayload();
    elements.assessButton.disabled = true;
    elements.assessButton.textContent = "Assessing…";
    const result = await runAssessment(payload);
    renderResult(result);
  } catch (error) {
    showError(shortError(error));
  } finally {
    elements.assessButton.disabled = false;
    elements.assessButton.textContent = "Assess event";
  }
});

elements.resetButton.addEventListener("click", () => {
  elements.form.reset();
  document.querySelectorAll("details").forEach((node) => {
    node.open = false;
  });
  clearError();
  clearResults();
});

elements.themeToggle.addEventListener("click", () => {
  const current = document.documentElement.dataset.theme;
  setTheme(current === "dark" ? "light" : "dark");
});

initializeTheme();
clearResults();
initializePython();
