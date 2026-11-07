"""Unit tests for the private formatter, independent of Formula parsing.

B and some isotope identities below are algorithmic fixtures, NOT
additions to molformula's supported-element or isotope registry.
"""

from itertools import permutations
from types import MappingProxyType

import pytest

from molformula._hill import _AtomKey, _format_hill


@pytest.mark.parametrize(
    "counts, expected",
    [
        ({("O", None): 1, ("H", None): 4, ("C", None): 1}, "CH4O"),
        ({("H", None): 8, ("O", None): 1, ("C", None): 3}, "C3H8O"),
        ({("N", None): 1, ("H", None): 7, ("C", None): 2}, "C2H7N"),
        ({("O", None): 6, ("C", None): 6, ("H", None): 12}, "C6H12O6"),
        ({("P", None): 4, ("O", None): 29, ("N", None): 15,
          ("H", None): 49, ("C", None): 38}, "C38H49N15O29P4"),
        ({("C", None): 1}, "C"),
        ({("H", None): 2, ("O", None): 1}, "H2O"),
        ({("N", None): 1, ("H", None): 3}, "H3N"),
        ({("S", None): 1, ("O", None): 4, ("H", None): 2}, "H2O4S"),
        # These cases distinguish true Hill order from a fixed C,H,... list.
        ({("H", None): 1, ("F", None): 1}, "FH"),
        ({("H", None): 3, ("F", None): 1, ("C", None): 1}, "CH3F"),
        ({("Cl", None): 1, ("C", None): 1, ("H", None): 3}, "CH3Cl"),
        ({("H", None): 1, ("Cl", None): 1}, "ClH"),
        ({("H", None): 3, ("B", None): 1}, "BH3"),
        ({("H", None): 3, ("B", None): 1, ("C", None): 1}, "CH3B"),
        ({("S", None): 1, ("P", None): 1, ("O", None): 1,
          ("N", None): 1, ("H", None): 1, ("C", None): 1}, "CHNOPS"),
    ],
    ids=[
        "methanol", "isopropanol-counts", "dimethylamine-counts", "glucose",
        "rna-example", "omit-count-one", "water", "no-carbon-ammonia",
        "no-carbon-H-O-S", "no-carbon-F-before-H", "carbon-H-before-F",
        "carbon-with-two-letter-symbol", "no-carbon-Cl-before-H",
        "no-carbon-B-before-H", "carbon-B-after-H", "all-ordinary-counts-one",
    ],
)
def test_hill_order_and_compact_format(counts, expected):
    assert _format_hill(counts) == expected


@pytest.mark.parametrize(
    "counts, expected",
    [
        ({("C", 13): 1, ("O", None): 1, ("H", None): 6,
          ("C", None): 1}, "C[13C]H6O"),
        ({("O", None): 6, ("H", None): 12, ("C", 13): 6}, "[13C]6H12O6"),
        ({("F", None): 1, ("H", None): 3, ("C", 13): 1}, "[13C]H3F"),
        ({("C", None): 0, ("F", None): 1, ("H", None): 3,
          ("C", 13): 1}, "[13C]H3F"),
        ({("H", 2): 2, ("O", None): 1}, "[2H]2O"),
        ({("H", 2): 1, ("O", None): 1, ("H", None): 1}, "H[2H]O"),
        ({("H", 2): 1, ("F", None): 1}, "F[2H]"),
        ({("H", 2): 1, ("C", None): 1, ("F", None): 1}, "C[2H]F"),
        ({("N", 15): 1, ("C", None): 2, ("H", None): 7}, "C2H7[15N]"),
        ({("O", 18): 1, ("O", None): 1, ("C", None): 1}, "CO[18O]"),
        ({("C", 13): 1}, "[13C]"),
        ({("C", 13): 2}, "[13C]2"),
        # Formatting-only fixtures: do NOT imply these nuclei are registered.
        ({("C", 13): 1, ("C", 12): 1, ("C", None): 1}, "C[12C][13C]"),
        ({("C", 13): 1, ("C", 9): 1, ("C", None): 1}, "C[9C][13C]"),
        ({("C", None): 1, ("O", 18): 1, ("N", 15): 1,
          ("N", None): 1, ("H", 2): 1, ("H", None): 1}, "CH[2H]N[15N][18O]"),
    ],
    ids=[
        "ordinary-C-before-C13", "explicit-carbon-only", "explicit-C-is-carbon",
        "zero-ordinary-C-with-positive-C13", "deuterated-water",
        "ordinary-H-before-deuterium", "no-carbon-with-deuterium",
        "carbon-with-deuterium", "nitrogen15", "ordinary-O-before-O18",
        "single-isotope", "isotope-count", "preserve-explicit-C12",
        "numeric-not-lexical-mass-order", "all-element-and-isotope-groups",
    ],
)
def test_isotope_rules(counts, expected):
    assert _format_hill(counts) == expected


@pytest.mark.parametrize(
    "counts, expected",
    [
        ({("C", None): 0, ("H", None): 1, ("F", None): 1}, "FH"),
        ({("C", 13): 0, ("H", None): 1, ("F", None): 1}, "FH"),
        ({("C", None): 1, ("H", None): 0, ("O", None): 1}, "CO"),
        ({("C", None): 1, ("C", 13): 0}, "C"),
    ],
)
def test_zero_entries_do_not_affect_order_or_output(counts, expected):
    assert _format_hill(counts) == expected


def test_insertion_order_does_not_change_output():
    items = [
        (("F", None), 1), (("H", 2), 1), (("C", 13), 1),
        (("H", None), 2), (("C", None), 1),
    ]
    for permuted in permutations(items):
        assert _format_hill(dict(permuted)) == "C[13C]H2[2H]F"


def test_does_not_mutate_input_and_accepts_read_only_mapping():
    counts: dict[_AtomKey, int] = {("O", None): 1, ("H", None): 4, ("C", None): 1, ("N", None): 0}
    before = list(counts.items())
    assert _format_hill(MappingProxyType(counts)) == "CH4O"
    assert list(counts.items()) == before


def test_large_integer_counts_are_not_rounded_or_abbreviated():
    assert _format_hill({("C", None): 12345678901234567890}) == "C12345678901234567890"


@pytest.mark.parametrize("counts", [{}, {("C", None): 0}, {("H", None): 0, ("O", 18): 0}])
def test_empty_formula_is_rejected(counts):
    with pytest.raises(ValueError, match="at least one atom"):
        _format_hill(counts)


@pytest.mark.parametrize("value", [-1, -100])
def test_negative_counts_are_never_silently_omitted(value):
    with pytest.raises(ValueError, match="non-negative"):
        _format_hill({("C", None): 1, ("H", None): value})


@pytest.mark.parametrize("value", [True, False, 1.0, "1", None])
def test_count_type_is_checked(value):
    with pytest.raises(TypeError, match="atom count"):
        _format_hill({("C", None): value})


@pytest.mark.parametrize("value", [True, 13.0, "13"])
def test_mass_number_type_is_checked(value):
    with pytest.raises(TypeError, match="mass number"):
        _format_hill({("C", value): 1})


@pytest.mark.parametrize("value", [0, -13])
def test_mass_number_must_be_positive(value):
    with pytest.raises(ValueError, match="mass number"):
        _format_hill({("C", value): 1})


@pytest.mark.parametrize("symbol", ["", "13C", "[13C]", "CH", "c", " C", "C\n"])
def test_only_resolved_element_symbol_shapes_are_accepted(symbol):
    with pytest.raises(ValueError, match="elemental symbol"):
        _format_hill({(symbol, None): 1})


@pytest.mark.parametrize("value", [None, "C", [(('C', None), 1)]])
def test_counts_must_be_a_mapping(value):
    with pytest.raises(TypeError, match="mapping"):
        _format_hill(value)


@pytest.mark.parametrize("atom", ["C", ("C",), ("C", None, 1)])
def test_key_shape_is_checked(atom):
    with pytest.raises(TypeError, match="atom keys"):
        _format_hill({atom: 1})


def test_symbol_type_is_checked():
    with pytest.raises(TypeError, match="elemental symbol"):
        _format_hill({(6, None): 1})  # type: ignore[arg-type]
