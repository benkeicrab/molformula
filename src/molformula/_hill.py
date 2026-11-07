# Copyright 2026 Hiroshi Nakayama and RIKEN
# SPDX-License-Identifier: BSD-3-Clause

"""Private, parser-independent Hill-style formula formatter.

Elemental order follows the Hill system.
    Edwin A. Hill, JACS (1900) 22 (8): 478-494.
    https://doi.org/10.1021/ja02046a005
Isotopes use this package's
additional rule: the ordinary symbol first, then explicit isotopes in
increasing mass-number order. No registry of supported atom species is
defined here. The caller supplies already resolved atom-species identities;
the parser and its Enum registry decide which atom species
the package supports.

This module is not a parser and must not be re-exported as public API.
It uses only Python 3.10-compatible syntax and standard-library features.
"""

from collections.abc import Mapping
import re
from typing import TypeAlias


# None means the ordinary symbol, NOT a mass number of zero.
_AtomKey: TypeAlias = tuple[str, int | None]
_ELEMENT_SYMBOL = re.compile(r"[A-Z][a-z]?")


def _hill_key(atom: _AtomKey, has_carbon: bool) -> tuple[int, str, int]:
    """Return an element-order key followed by an isotope-order key."""
    symbol, mass_number = atom
    if has_carbon:
        priority = 0 if symbol == "C" else 1 if symbol == "H" else 2
    else:
        priority = 0
    return (priority, symbol, 0 if mass_number is None else mass_number)


def _format_hill(counts: Mapping[_AtomKey, int]) -> str:
    """Format resolved atom counts as a compact formula string.

    Counts of one are omitted. Explicit isotopes are written in bracket
    notation, e.g. ``[13C]``. An explicit carbon isotope also makes
    a formula carbon-containing for Hill ordering.

    Keys are (elemental_symbol, mass_number_or_None), e.g. ("C", None)
    for ordinary C, and ("C", 13) for explicitly labeled 13C.
    Deuterium MUST be supplied as ("H", 2), not as the independent
    elemental symbol "D".

    Input may include zero-count entries, which are omitted. Negative
    counts, non-integer counts (including bool), and an empty/all-zero
    formula are rejected. Input is never modified. Explicit isotope
    labels are retained rather than merged with an ordinary symbol.

    This function checks internal data shape, not Enum membership,
    atomic masses, isotope stability, chemical valence, or molecular
    structure. Invalid internal input raises built-in TypeError
    /ValueError, not public parser or formula-operation exceptions.
    """
    if not isinstance(counts, Mapping):
        raise TypeError("counts must be a mapping")

    positive: dict[_AtomKey, int] = {}
    for atom, count in counts.items():
        if type(atom) is not tuple or len(atom) != 2:
            raise TypeError("atom keys must be (elemental_symbol, mass_number) tuples")
        symbol, mass_number = atom
        if type(symbol) is not str:
            raise TypeError("elemental symbol must be a string")
        if _ELEMENT_SYMBOL.fullmatch(symbol) is None:
            raise ValueError("invalid elemental symbol shape")
        if mass_number is not None:
            if type(mass_number) is not int:
                raise TypeError("mass number must be an integer or None")
            if mass_number <= 0:
                raise ValueError("mass number must be positive")
        if type(count) is not int:
            raise TypeError("atom count must be an integer")
        if count < 0:
            raise ValueError("atom count must be non-negative")
        if count:
            positive[atom] = count

    if not positive:
        raise ValueError("formula must contain at least one atom")

    # A positive explicit-carbon-isotope count also establishes carbon.
    has_carbon = any(symbol == "C" for symbol, _ in positive)
    ordered_atoms = sorted(positive,
                           key=lambda atom: _hill_key(atom, has_carbon))

    parts: list[str] = []
    for symbol, mass_number in ordered_atoms:
        atom = (symbol, mass_number)
        token = symbol if mass_number is None else f"[{mass_number}{symbol}]"
        count = positive[atom]
        parts.append(token if count == 1 else f"{token}{count}")
    return "".join(parts)
