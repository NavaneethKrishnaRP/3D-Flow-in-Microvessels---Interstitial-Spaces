# Contributing

Thanks for considering a contribution!

## Dev setup
- MATLAB (R2022b+ recommended), with required toolboxes if applicable.
- Java 21 + Maven (or switch to Gradle if you prefer).

## Running tests
- MATLAB: `matlab -batch "runtests"` (tests live in `matlab/tests`).
- Java: `mvn -q verify` (JUnit tests live in `java/src/test/java`).

## Branching & commits
- Use feature branches like `feat/...`, `fix/...`, `docs/...`.
- Prefer Conventional Commits style (e.g., `feat: add FEM postprocessing`).
- Open a PR; ensure CI passes and update docs if behavior changes.

## Code style
- MATLAB: functions under `matlab/functions`, scripts under `matlab/scripts`.
- Java: follow standard package + class structure in `java/src/main/java`.

## Data
- Commit only *small* sample data in `data/sample`.
- Large files should be stored externally; document how to fetch them in `data/README.md`.
