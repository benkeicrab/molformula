# Copyright 2026 Hiroshi Nakayama and RIKEN
# SPDX-License-Identifier: BSD-3-Clause

"""Parser for molecular formulas and simple condensed formulas.

Grammar summary:
- Elemental symbols consist of one uppercase ASCII letter optionally
  followed by one lowercase ASCII letter.
- Explicit isotopes are written as [13C], [2H], etc.
- Atom-species counts are positive integers; 1 may be omitted.
- Repeated atom species are summed, e. g. CH3OH -> H: 3 + 1 = 4
- Parentheses may be nested and may be followed immediately by a positive multiplier.
- Whitespace may separate tokens but may not split lexical tokens.
- Only atom species contained in ``supported_atom_species`` are accepted.
- Complexes, component separators, charges, and structural notation
  are not supported, e. g. K3[Fe(CN)6], CuSO4.H2O.

The parser returns atom-species counts only. Input ordering and
grouping are not preserved.
"""

from dataclasses import dataclass
from typing import AbstractSet
from collections import Counter

from . import _error


@dataclass(slots=True)
class _Group:
    counts: Counter[str]
    opening_position: int


def _skip_whitespace(text: str, position: int) -> int:
    while position < len(text) and text[position].isspace():
        position += 1
    return position


def _match_symbol(text: str, position: int) -> tuple[str | None, int]:
    start = _skip_whitespace(text, position)

    if start >= len(text) or text[start] not in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
        return None, start

    candidate = text[start]
    next_pos = start + 1

    if next_pos < len(text) and text[next_pos] in 'abcdefghijklmnopqrstuvwxyz':
        candidate += text[next_pos]
        return candidate, next_pos + 1

    return candidate, start + 1


def _parse_element_symbol(text: str, position: int) -> tuple[str, int]:
    candidate, next_position = _match_symbol(text, position)
    if candidate is None:
        raise _error.InvalidFormulaError(next_position)

    return candidate, next_position


def _require_supported(
        supported_atom_species: AbstractSet[str],
        symbol: str, position: int
        ) -> None:
    if symbol not in supported_atom_species:
        raise _error.UnsupportedAtomSpeciesError(symbol, position)


def _match_int(text: str, position: int) -> tuple[int | None, int]:
    digits: list[str] = []

    while True:
        if position >= len(text) or text[position] not in "0123456789":
            break

        digits.append(text[position])
        position += 1

    if not digits:
        return None, position

    return int("".join(digits)), position


def _parse_count(text: str, position: int) -> tuple[int, int]:
    start = _skip_whitespace(text, position)

    if start < len(text) and text[start] in '-0':
        raise _error.InvalidCountError(text[start], start)

    count, next_position = _match_int(text, start)
    if count is None:
        count = 1

    return count, next_position


def _parse_isotope_symbol(
        text: str,
        position: int
        ) -> tuple[str, int]:
    start = _skip_whitespace(text, position)

    if start >= len(text) or text[start] not in '123456789':
        raise _error.InvalidBracketError(start)

    mass_num, next_position = _match_int(text, start)

    symbol_start = next_position

    elem_symbol, next_position = _match_symbol(text, symbol_start)
    if elem_symbol is None:
        raise _error.InvalidBracketError(symbol_start)

    isotope_symbol = str(mass_num) + elem_symbol

    return isotope_symbol, next_position


def parse_formula_str(
        supported_atom_species: AbstractSet[str],
        formula_str: str
        ) -> dict[str, int]:
    """Parse a formula string and return aggregated atom-species counts.

    `atom_species` contains the supported atom-species symbols, using
    unbracketed isotope names such as `"13C"` and `"2H"`.

    Raises:
        TypeError: If `formula_str` is not a string.
        EmptyParenthesisError: If a parenthesized group is empty.
        UnmatchedParenthesesError: If an opening or closing parenthesis is unmatched.
        InvalidBracketError: If bracketed isotope notation is malformed,
            including a missing or unexpected bracket, mass number, or
            element symbol.
        EmptyFormulaError: If the input is empty or contains only whitespace.
        InvalidFormulaError: If an atom species cannot be parsed.
        UnsupportedAtomSpeciesError: If an atom species is not among the supported species.
        InvalidCountError: If an atom-species count or group multiplier is zero or negative.

    Examples:
    >>> parse_formula_str(frozenset({'C', '13C', 'H', 'O'}), 'C5[13C]H12O6')
    {'C': 5, '13C': 1, 'H': 12, 'O': 6}
    >>> parse_formula_str(frozenset({'C', '13C', 'H', 'O'}), '  ')
    Traceback (most recent call last):
        ...
    molformula._error.EmptyFormulaError: Empty formula
    """
    # if not isinstance(formula_str, str):
    if type(formula_str) is not str:  # more strict
        raise TypeError("formula must be a string")

    if _skip_whitespace(formula_str, 0) >= len(formula_str):
        raise _error.EmptyFormulaError(0)

    i = 0
    stack = [_Group(Counter(), -1)]
    length = len(formula_str)

    while i < length:
        position = _skip_whitespace(formula_str, i)
        i = position
        if i >= length:
            break

        c = formula_str[i]

        if c == '(':
            opening_position = i
            stack.append(_Group(Counter(), opening_position))
            i += 1
            continue

        elif c == ')':
            closing_position = i
            i += 1

            if len(stack) == 1:
                raise _error.UnmatchedParenthesesError(closing_position)

            top = stack.pop()
            if not top.counts:
                raise _error.EmptyParenthesisError(top.opening_position)

            multiplier, next_position = _parse_count(formula_str, i)
            i = next_position

            curr = stack[-1]
            for elm, cnt in top.counts.items():
                curr.counts[elm] += cnt * multiplier
            continue

        elif c == '[':
            bracket_start = i
            candidate, next_position = _parse_isotope_symbol(formula_str, bracket_start + 1)
            i = _skip_whitespace(formula_str, next_position)

            if i >= length or formula_str[i] != ']':
                # Report exact invalid position.
                raise _error.InvalidBracketError(i)
            i += 1

            isotope_start = _skip_whitespace(formula_str, bracket_start + 1)
            _require_supported(supported_atom_species, candidate, isotope_start)
            elem = candidate

        elif c == ']':
            raise _error.InvalidBracketError(i)

        else:
            symbol_start = i
            candidate, next_position = _parse_element_symbol(formula_str, symbol_start)
            _require_supported(supported_atom_species, candidate, symbol_start)
            elem = candidate
            i = next_position

        count, next_position = _parse_count(formula_str, i)
        i = next_position

        curr = stack[-1]
        curr.counts[elem] += count

    if len(stack) > 1:
        raise _error.UnmatchedParenthesesError(
            stack[-1].opening_position
        )

    return dict(stack[0].counts)
