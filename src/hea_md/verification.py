"""Recalculate the campaign and compare all exported numerical representations."""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

from .campaign import calculate, resolve_inside
from .io import sha256_file

EXPORTS = {
    'summary.csv': 'summary',
    'sensitivity.csv': 'sensitivity',
    'method-sensitivity.csv': 'method_sensitivity_ranges',
    'paired-sensitivity.csv': 'paired_sensitivity',
}


def compare_results(expected: Any, actual: Any, path: str = 'root') -> None:
    """Compare structures exactly and finite floats within calculation tolerance."""
    if isinstance(expected, dict):
        if not isinstance(actual, dict) or set(expected) != set(actual):
            raise ValueError(f'{path}: keys differ.')
        for key in expected:
            compare_results(expected[key], actual[key], f'{path}.{key}')
    elif isinstance(expected, list):
        if not isinstance(actual, list) or len(expected) != len(actual):
            raise ValueError(f'{path}: length differs.')
        for i, (left, right) in enumerate(zip(expected, actual)):
            compare_results(left, right, f'{path}[{i}]')
    elif isinstance(expected, float):
        if (isinstance(actual, bool) or not isinstance(actual, (int, float))
                or not math.isfinite(expected) or not math.isfinite(actual)
                or not math.isclose(expected, actual, rel_tol=1e-10, abs_tol=1e-11)):
            raise ValueError(f'{path}: numerical result differs.')
    elif type(expected) is not type(actual) or expected != actual:
        raise ValueError(f'{path}: values differ.')


def compare_csv(path: Path, expected: list[dict]) -> None:
    """Check headers, ordering, row counts, numerical cells and empty fields."""
    with path.open(encoding='utf-8', newline='') as stream:
        reader = csv.reader(stream)
        header = next(reader, None)
        if header != list(expected[0]):
            raise ValueError(f'{path.name}: column names or ordering differ.')
        rows = list(reader)
    if len(rows) != len(expected):
        raise ValueError(f'{path.name}: row count differs.')
    for index, (row, reference) in enumerate(zip(rows, expected), 2):
        if len(row) != len(header):
            raise ValueError(f'{path.name}:{index}: column count differs.')
        for field, text in zip(header, row):
            value = reference[field]
            location = f'{path.name}:{index}:{field}'
            if value is None:
                if text != '':
                    raise ValueError(f'{location}: expected an empty field.')
            elif isinstance(value, float):
                try:
                    number = float(text)
                except ValueError as exc:
                    raise ValueError(f'{location}: invalid number.') from exc
                compare_results(value, number, location)
            elif text != str(value):
                raise ValueError(f'{location}: value differs.')


def verify(root: Path, results: Path | None = None) -> dict:
    """Check source-data hashes, full JSON, and four CSV exports."""
    root = Path(root).resolve()
    results = Path(results) if results is not None else root / 'results/analysis.json'
    if not results.is_absolute():
        results = root / results
    computed = calculate(root)
    stored = json.loads(results.read_text(encoding='utf-8'))
    compare_results(computed, stored)
    for filename, key in EXPORTS.items():
        compare_csv(results.parent / filename, computed[key])
    return computed
