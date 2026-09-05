"""Parsing helpers shared by the tutorial regression tests."""
import csv
import re

NUM_RE = re.compile(r"[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?")


def extract_numbers(text):
    """Pull every float literal out of a string (works whether it's '[ 1 2 3 ]' or '1 2 3')."""
    return [float(x) for x in NUM_RE.findall(text)]


def parse_comparison_pairs(stdout):
    """Pair up ChemGen/Cantera comparison lines emitted by default_test.py.

    default_test.py always prints a "Chemgen ... : <values>" (or "Source test
    result:" / "ChemGen internal energy:") line immediately followed by the
    matching "Cantera ... : <values>" line. Returns a list of
    (label, chemgen_values, cantera_values).
    """
    pairs = []
    pending_label = None
    pending_values = None
    for line in stdout.splitlines():
        stripped = line.strip()
        if not stripped or ":" not in stripped:
            continue
        label, _, rest = stripped.partition(":")
        label = label.strip()
        if label.startswith(("Chemgen", "ChemGen", "Source test result")):
            pending_label = label
            pending_values = extract_numbers(rest)
        elif label.startswith("Cantera") and pending_values is not None:
            pairs.append((pending_label, pending_values, extract_numbers(rest)))
            pending_label = None
            pending_values = None
    return pairs


def max_relative_error(chemgen_values, cantera_values, floor=1e-8):
    """Relative error, falling back to absolute error for near-zero references."""
    assert len(chemgen_values) == len(cantera_values), (
        f"length mismatch: chemgen has {len(chemgen_values)}, cantera has {len(cantera_values)}"
    )
    errors = []
    for cg, ct_val in zip(chemgen_values, cantera_values):
        denom = abs(ct_val) if abs(ct_val) > floor else 1.0
        errors.append(abs(cg - ct_val) / denom)
    return max(errors) if errors else 0.0


def read_l2_norm_csv(path):
    """Read tutorial 02's l2_norm_results.csv -> list of (temperature, l2_norm) floats."""
    rows = []
    with open(path, newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        assert len(header) == 2, f"unexpected l2_norm_results.csv header: {header}"
        for row in reader:
            if not row:
                continue
            rows.append((float(row[0]), float(row[1])))
    return rows
