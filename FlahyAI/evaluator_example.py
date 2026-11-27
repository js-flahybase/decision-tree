# This is an example of structure of an evaluator function along with its schema.
# All connected dependencies and their working are available in *_AI.py files.

def evaluate_type2_diabetes(labs, patient, symptoms, family_history, past_history):
    # 2. Pull values from labs / patient
    fpg = labs.get("fasting_glucose")
    hba1c = labs.get("hba1c")
    # eag = labs.get("estimated_average_glucose_(eag)")
    fasting_insulin = labs.get("fasting_insulin")

    if fpg is None or hba1c is None or fasting_insulin is None:
        return [{"Condition": "Type 2 Diabetes", "Category": INSUFFICIENT_DATA, "Values": "", "TriggeredValues": "", "PartialTriggered": ""}]

    # 3. Check if values are abnormal
    fpg_flag = is_elevated(fpg, THRESHOLDS["fasting_plasma_glucose_high"])
    hba1c_flag = is_elevated(hba1c, THRESHOLDS["hba1c_high"])
    # eag_flag = is_elevated(eag, THRESHOLDS["eag_high"])
    insulin_flag = is_above(fasting_insulin, THRESHOLDS["fasting_insulin_high"])

    fpg_borderline_flag = in_range(fpg, THRESHOLDS2["fasting_plasma_glucose_borderline_low"], THRESHOLDS2["fasting_plasma_glucose_borderline_high"])
    hba1c_borderline_flag = in_range(hba1c, THRESHOLDS2["hba1c_borderline_low"], THRESHOLDS2["hba1c_borderline_high"])
    insulin_borderline_flag=is_above(fasting_insulin,THRESHOLDS2["fasting_insulin_borderline_high"])
    # 4. Decide category
    triggered = []
    partial_triggered = []
    if fpg_flag and hba1c_flag and insulin_flag or (fpg_flag and hba1c_flag and insulin_flag and (symptoms or family_history or past_history)):
        category = "Significant Pattern"
        _track(triggered, "fasting_glucose", fpg, THRESHOLDS["fasting_plasma_glucose_high"], ">=", fpg_flag)
        _track(triggered, "hba1c", hba1c, THRESHOLDS["hba1c_high"], ">=", hba1c_flag)
        _track(triggered, "fasting_insulin", fasting_insulin, THRESHOLDS["fasting_insulin_high"], ">", insulin_flag)
    elif fpg_borderline_flag and hba1c_borderline_flag and insulin_borderline_flag or (fpg_borderline_flag and hba1c_borderline_flag and insulin_borderline_flag and (symptoms or family_history or past_history)):
        category = "Early Pattern"
        _track_range(triggered, "fasting_glucose", fpg, THRESHOLDS2["fasting_plasma_glucose_borderline_low"], THRESHOLDS2["fasting_plasma_glucose_borderline_high"], fpg_borderline_flag)
        _track_range(triggered, "hba1c", hba1c, THRESHOLDS2["hba1c_borderline_low"], THRESHOLDS2["hba1c_borderline_high"], hba1c_borderline_flag)
        _track(triggered, "fasting_insulin", fasting_insulin, THRESHOLDS2["fasting_insulin_borderline_high"], ">", insulin_borderline_flag)
    else:
        category = "Typical"
        if fpg_flag:
            _track(partial_triggered, "fasting_glucose", fpg, THRESHOLDS["fasting_plasma_glucose_high"], ">=", fpg_flag)
        elif fpg_borderline_flag:
           _track_range(partial_triggered, "fasting_glucose", fpg, THRESHOLDS2["fasting_plasma_glucose_borderline_low"], THRESHOLDS2["fasting_plasma_glucose_borderline_high"], fpg_borderline_flag)
        if hba1c_flag:
            _track(partial_triggered, "hba1c", hba1c, THRESHOLDS["hba1c_high"], ">=", hba1c_flag)
        elif hba1c_borderline_flag:
            _track_range(partial_triggered, "hba1c", hba1c, THRESHOLDS2["hba1c_borderline_low"], THRESHOLDS2["hba1c_borderline_high"], hba1c_borderline_flag)
        # _track(partial_triggered, "eag", eag, THRESHOLDS["eag_high"], ">=", eag_flag)
        if insulin_flag:
            _track(partial_triggered, "fasting_insulin", fasting_insulin, THRESHOLDS["fasting_insulin_high"], ">", insulin_flag)
        elif insulin_borderline_flag:
            _track(partial_triggered, "fasting_insulin", fasting_insulin, THRESHOLDS2["fasting_insulin_borderline_high"], ">", insulin_borderline_flag)

    return [{
        "Condition": "Type 2 Diabetes",
        "Category": category,
        "Values": _fmt_values(fasting_glucose=fpg, hba1c=hba1c, fasting_insulin=fasting_insulin),
        "TriggeredValues": ", ".join(triggered),
        "PartialTriggered": ", ".join(partial_triggered)
    }]
evaluate_type2_diabetes_schema = {
    'type': 'function',
    'function': {
        'name': 'evaluate_type2_diabetes',
        'description': 'Evaluate Type 2 Diabetes susceptibility from fasting glucose, hba1c and fasting insulin. Blood-only - no genetics dependency.',
        'parameters': {
            'type': 'object',
            'properties': {
                'fasting_glucose': {
                    'type': ['number', 'null'],
                    'description': 'fasting_glucose, mg/dL. Null if not tested.',
                },
                'hba1c': {
                    'type': ['number', 'null'],
                    'description': 'hba1c, %. Null if not tested.',
                },
                'fasting_insulin': {
                    'type': ['number', 'null'],
                    'description': 'fasting_insulin, uIU/mL. Null if not tested.',
                },
                'family_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Type 2 Diabetes', {}).get('family_history', []),
                    },
                    'description': 'Which of these the patient reports in family history.',
                },
                'past_history': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Type 2 Diabetes', {}).get('past_history', []),
                    },
                    'description': "Which of these are in the patient's own past medical history.",
                },
                'symptoms': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': USER_CONTEXT.get('Type 2 Diabetes', {}).get('symptoms', []),
                    },
                    'description': 'Which of these symptoms the patient reports.',
                },
            },
            'required': ['fasting_glucose', 'hba1c', 'fasting_insulin', 'family_history', 'past_history', 'symptoms'],
            'additionalProperties': False,
        },
        'strict': True,
    },
}