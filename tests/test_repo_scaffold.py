from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).parent.parent


def test_license_first_line_is_mit():
    first_line = (REPO_ROOT / "LICENSE").read_text().splitlines()[0]
    assert first_line == "MIT License"


def test_ci_workflow_has_lint_and_test_steps():
    workflow = yaml.safe_load((REPO_ROOT / ".github" / "workflows" / "ci.yml").read_text())
    jobs = workflow["jobs"]
    steps = [step for job in jobs.values() for step in job["steps"]]
    run_commands = " ".join(step.get("run", "") for step in steps)
    assert "ruff" in run_commands
    assert "pytest" in run_commands
