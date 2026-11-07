"""Acceptance tests for the formula grammar.
"""
from typing import AbstractSet, Any

import pytest

from molformula._parser import parse_formula_str
import molformula._error as _error
from molformula._core import _INDEX_BY_SYMBOL, _rep_from_counts


VALID_CASES = (
    pytest.param("C6H12O6", {"C": 6, "H": 12, "O": 6}, id="molecular-formula"),
    pytest.param("CH3OH", {"C": 1, "H": 4, "O": 1}, id="repeated-hydrogen"),
    pytest.param("C5[13C]H12O6", {"C": 5, "13C": 1, "H": 12, "O": 6}, id="carbon-isotope"),
    pytest.param("CH[2H]2OH", {"C": 1, "H": 2, "2H": 2, "O": 1}, id="hydrogen-isotope"),
    pytest.param("(CH3)2CHOH", {"C": 3, "H": 8, "O": 1}, id="group-multiplier"),
    pytest.param("( CH3 ) 2CHOH", {"C": 3, "H": 8, "O": 1}, id="whitespace-in-group"),
    pytest.param("CH3CH(OH)CH3", {"C": 3, "H": 8, "O": 1}, id="group-default-one"),
    pytest.param("((CH3)2CH)2O", {"C": 6, "H": 14, "O": 1}, id="nested-groups"),
    pytest.param("(C6H12O6)", {"C": 6, "H": 12, "O": 6}, id="outer-group"),
    pytest.param("C6 H12 O6", {"C": 6, "H": 12, "O": 6}, id="whitespace-between-terms"),
    pytest.param("[ 13 C ]6H12O6", {"13C": 6, "H": 12, "O": 6}, id="whitespace-in-isotope1"),
    pytest.param("CCC", {"C": 3}, id="repeated-element"),
    pytest.param("([13C]H3)2NH", {"13C": 2, "H": 7, "N": 1}, id="combination-rules")
)

INVALID_CASES = (
    pytest.param("-H2O", _error.InvalidFormulaError, id="leading-minus"),
    pytest.param("2H2O", _error.InvalidFormulaError, id="leading-number"),
    pytest.param(" 2H2O", _error.InvalidFormulaError, id="leading-space-number"),
    pytest.param("CH3@CC", _error.InvalidFormulaError, id="illegal-character"),
    pytest.param("H2O.H2O", _error.InvalidFormulaError, id="illegal-character2"),
    pytest.param("H2O·H2O", _error.InvalidFormulaError, id="illegal-character3"),
    pytest.param("C6 H1 2O6", _error.InvalidFormulaError, id="whitespace-inside-count"),
    pytest.param("NH4+", _error.InvalidFormulaError, id="illegal-character4"),
    pytest.param("CH3-CH2-OH", _error.InvalidFormulaError, id="illegal-character5"),
    pytest.param("", _error.EmptyFormulaError, id="empty-input"),
    pytest.param("   \t\n", _error.EmptyFormulaError, id="tab-space-eol-input"),
    pytest.param("(CH32CHOH", _error.UnmatchedParenthesesError, id="unclosed-group"),
    pytest.param("CH3OH)", _error.UnmatchedParenthesesError, id="unopened-group"),
    pytest.param(")CH32CHOH", _error.UnmatchedParenthesesError, id="unopened-group2"),
    pytest.param(")0", _error.UnmatchedParenthesesError, id="unopened-group3"),
    pytest.param("((H)2))", _error.UnmatchedParenthesesError, id="unopened-nested-group"),
    pytest.param("()CH32CHOH", _error.EmptyParenthesisError, id="empty-group1"),
    pytest.param("()0", _error.EmptyParenthesisError, id="empty-group2"),
    pytest.param("GaOH", _error.UnsupportedAtomSpeciesError, id="unregistered-element"),
    pytest.param("N aOH", _error.InvalidFormulaError, id="unregistered-element"),
    pytest.param("[12C]6H12O6", _error.UnsupportedAtomSpeciesError, id="unregistered-isotope"),
    pytest.param("C06H12O6", _error.InvalidCountError, id="leading-zero-count1"),
    pytest.param("(CH3)02", _error.InvalidCountError, id="leading-zero-count2"),
    pytest.param("C 06H12O6", _error.InvalidCountError, id="leading-zero-count1"),
    pytest.param("(CH3) 02", _error.InvalidCountError, id="leading-zero-count2"),
    pytest.param("C0", _error.InvalidCountError, id="zero-count1"),
    pytest.param("(CH3)0", _error.InvalidCountError, id="zero-count2"),
    pytest.param("[13C]0", _error.InvalidCountError, id="zero-count3"),
    pytest.param("C 0", _error.InvalidCountError, id="zero-count1"),
    pytest.param("(CH3) 0", _error.InvalidCountError, id="zero-count2"),
    pytest.param("[13C] 0", _error.InvalidCountError, id="zero-count3"),
    pytest.param("C-6H12O6", _error.InvalidCountError, id="negative-count1"),
    pytest.param("C -6H12O6", _error.InvalidCountError, id="negative-count1"),
    pytest.param("(CH3)-2", _error.InvalidCountError, id="negative-count2"),
    pytest.param("[13C]-2", _error.InvalidCountError, id="negative-count3"),
    pytest.param("(CH3) -2", _error.InvalidCountError, id="negative-count2"),
    pytest.param("[13C] -2", _error.InvalidCountError, id="negative-count3"),
    pytest.param("[CH3OH]", _error.InvalidBracketError, id="non-isotope-bracket"),
    pytest.param("[0C]6H12O6", _error.InvalidBracketError, id="zero-mass-number"),
    pytest.param("[013C]6H12O6", _error.InvalidBracketError, id="leading-zero-mass-number"),
    pytest.param("[ 1 3C ]6H12O6", _error.InvalidBracketError, id="whitespace-in-isotope-symbol"),
    pytest.param("[[13C]]6H12O6", _error.InvalidBracketError, id="nested-bracket1"),
    pytest.param("[13[C]]6H12O6", _error.InvalidBracketError, id="nested-bracket2"),
    pytest.param("[13C]]6H12O6", _error.InvalidBracketError, id="unmatched-bracket1"),
    pytest.param("[", _error.InvalidBracketError, id="unmatched-bracket2"),
    pytest.param("[13", _error.InvalidBracketError, id="unmatched-bracket3"),
    pytest.param("[13C", _error.InvalidBracketError, id="unmatched-bracket4"),
    pytest.param("[]", _error.InvalidBracketError, id="empty-isotope"),
    pytest.param("[13]", _error.InvalidBracketError, id="incomplete-isotope"),
    pytest.param("[13c]", _error.InvalidBracketError, id="invalid-isotope1"),
    pytest.param("[13@]", _error.InvalidBracketError, id="invalid-isotope2"),
    pytest.param("[12C6H12O6", _error.InvalidBracketError, id="unregistered-isotope"),
)

TWO_LETTER_VALID_CASES = (
    pytest.param(
        frozenset({"Na", "H", "O"}),
        "NaOH", {"Na": 1, "H": 1, "O": 1},
        id="Na-supported",
    ),
    pytest.param(
        frozenset({"23Na", "H", "O"}),
        "[23Na]OH", {"23Na": 1, "H": 1, "O": 1},
        id="23Na-supported",
    ),
    pytest.param(
        frozenset({"23Na", "H", "O"}),
        "[23 Na]OH", {"23Na": 1, "H": 1, "O": 1},
        id="23 Na-supported",
    ),
    pytest.param(
        frozenset({"Co", "C", "O"}),
        "Co", {"Co": 1},
        id="Co-supported",
    ),
    pytest.param(
        frozenset({"Co", "C", "O"}),
        "CO", {"C": 1, "O": 1},
        id="CO-is-not-Co",
    ),
)

TWO_LETTER_INVALID_CASES = (
    pytest.param(
        frozenset({"N", "H", "O"}),
        "NaOH", _error.UnsupportedAtomSpeciesError,
        id="Na-unsupported",
    ),
    pytest.param(
        frozenset({"Na", "H", "O"}),
        "[23Na]OH", _error.UnsupportedAtomSpeciesError,
        id="23Na-unsupported",
    ),
    pytest.param(
        frozenset({"23Na", "H", "O"}),
        "[2 3Na]OH", _error.InvalidBracketError,
        id="2 3Na-unsupported",
    ),
    pytest.param(
        frozenset({"23Na", "H", "O"}),
        "[23N a]OH", _error.InvalidBracketError,
        id="23N a-unsupported",
    ),
)

POSITION_CASES = (
    ("C6 H12 O-1", _error.InvalidCountError, 8),
    ("C6 @H12", _error.InvalidFormulaError, 3),
    ("C6 )", _error.UnmatchedParenthesesError, 3),
    ("( C6", _error.UnmatchedParenthesesError, 0),
    ("C[ 12C]H", _error.UnsupportedAtomSpeciesError, 3),
    ("C[13 C", _error.InvalidBracketError, 6),
)

_ELEMENTS = frozenset({'C', '13C', 'H', '2H', 'N', 'O'})


@pytest.mark.parametrize("text, expected", VALID_CASES)
def test_parser_valid_examples(
        text: str, expected: dict[str, int]
    ) -> None:
    formula_dict = parse_formula_str(_ELEMENTS, text)
    assert formula_dict == expected


@pytest.mark.parametrize("text, expected_error", INVALID_CASES)
def test_parser_invalid_examples(
        text: str, expected_error: Any
    ) -> None:
    with pytest.raises(expected_error) as caught:
        parse_formula_str(_ELEMENTS, text)
    assert caught.type is expected_error


@pytest.mark.parametrize("elements, text, expected", TWO_LETTER_VALID_CASES)
def test_element_symbol_parsing_valid(
        elements: AbstractSet[str], text: str, expected: dict[str, int]
    ) -> None:
    formula_dict = parse_formula_str(elements, text)
    assert type(formula_dict) is dict
    assert formula_dict == expected


@pytest.mark.parametrize("elements, text, expected_error", TWO_LETTER_INVALID_CASES)
def test_element_symbol_parsing_invalid(
        elements: AbstractSet[str], text: str, expected_error: Any,
    ) -> None:
    with pytest.raises(expected_error) as caught:
        parse_formula_str(elements, text)
    assert caught.type is expected_error


def test_non_str_input() -> None:
    with pytest.raises(TypeError):
        parse_formula_str(
            _ELEMENTS,
            ["H", "2", "O"]  # pyright: ignore[reportArgumentType]
        )


def test_unknown_element_error_attributes() -> None:
    with pytest.raises(_error.UnsupportedAtomSpeciesError) as caught:
        parse_formula_str(_ELEMENTS, "C NaH")

    assert caught.value.symbol == "Na"
    assert caught.value.position == 2


def test_unknown_isotope_error_attributes() -> None:
    with pytest.raises(_error.UnsupportedAtomSpeciesError) as caught:
        parse_formula_str(_ELEMENTS, "C[12C]H")

    assert caught.value.symbol == "12C"
    assert caught.value.position == 2


def test_rep_from_counts_with_isotope() -> None:
    counts = {
        "C": 5,
        "13C": 1,
        "H": 12,
        "O": 6,
    }

    rep = _rep_from_counts(counts)

    assert rep._counts[_INDEX_BY_SYMBOL["C"]] == 5
    assert rep._counts[_INDEX_BY_SYMBOL["13C"]] == 1
    assert rep._counts[_INDEX_BY_SYMBOL["H"]] == 12
    assert rep._counts[_INDEX_BY_SYMBOL["O"]] == 6
    assert rep._counts[_INDEX_BY_SYMBOL["N"]] == 0


def test_parser_result_converts_to_rep() -> None:
    counts = parse_formula_str(
        _ELEMENTS,
        "C5[13C]H12O6",
    )

    rep = _rep_from_counts(counts)

    assert rep._counts[_INDEX_BY_SYMBOL["C"]] == 5
    assert rep._counts[_INDEX_BY_SYMBOL["13C"]] == 1
    assert rep._counts[_INDEX_BY_SYMBOL["H"]] == 12
    assert rep._counts[_INDEX_BY_SYMBOL["O"]] == 6
    assert rep._counts[_INDEX_BY_SYMBOL["N"]] == 0


@pytest.mark.parametrize("text, expected_error, position", POSITION_CASES)
def test_invalid_positions(
    text: str, expected_error: Any, position: int,
) -> None:
    with pytest.raises(expected_error) as caught:
        parse_formula_str(_ELEMENTS, text)
    assert caught.type is expected_error
    assert caught.value.position == position
