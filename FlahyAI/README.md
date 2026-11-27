# FlahyAI — Structure Overview

This folder is a trimmed-down snapshot of how FlahyAI evaluates a patient's
data and talks to them about it.

Two layers are involved:

- **The evaluators** (Python) — pure decision-tree logic. Given a patient's
  recorded data, they return a diagnosis/category and recommendations. They
  don't talk to the model or the patient directly.
- **The agent/prompt layer** — the instructions that tell the LLM *how* to
  behave: when to ask questions, when to call an evaluator, how to handle
  corrections, how to phrase results. This is split across a system prompt
  and the tool-calling setup at the API level.

## Files in this folder

- **`agent.ts`** — the developer-facing API calling layer: where the technical-level system prompt lives. 

- **`Prompt_extra.txt`** — What sequence to follow before reaching an evaluator and how to perform interpretations.

- **`evaluator_example.py`** — It is a minimal, self-contained example (a  blood-based T2D evaluator) showing the two halves of one evaluator structure:
  1. a **schema** describing the inputs it expects (lab values, patient fields) and their types/units, and
  2. the **function** that takes those inputs and returns a decision.

- **`blood_conditions_functions_AI.py` and `blood_conditions_functions_SDT_AI.py`** — the full evaluator functions for
  all conditions, written for the live AI-driven flow (Just for reference).

