# Copyright 2026 Hiroshi Nakayama and RIKEN
# SPDX-License-Identifier: BSD-3-Clause

from importlib.metadata import version

from ._core import (
    Formula, IonFormula,
    atom_species_counts, exact_mass,
    protonated, deprotonated, deelectronated, ion_mz, apply_formula_change,
    make_formula_change, sum_formulas, to_str,
)
from ._error import (
    ChemicalFormulaError, FormulaChangeError, InvalidProtonCountError,
    IonizationError, InsufficientHydrogenError, EmptyIonFormulaError,
    FormulaParseError, EmptyParenthesisError, UnmatchedParenthesesError,
    InvalidBracketError, EmptyFormulaError, InvalidFormulaError,
    UnsupportedAtomSpeciesError, InvalidCountError,
)

__version__ = version(__package__)  # type: ignore[arg-type]

__all__ = [
    'Formula',
    'IonFormula',
    'to_str',
    'atom_species_counts',
    'exact_mass',
    'protonated',
    'deprotonated',
    'deelectronated',
    'ion_mz',
    'apply_formula_change',
    'make_formula_change',
    'sum_formulas',
    'ChemicalFormulaError',
    'FormulaChangeError',
    'InvalidProtonCountError',
    'IonizationError',
    'InsufficientHydrogenError',
    'EmptyIonFormulaError',
    'FormulaParseError',
    'EmptyParenthesisError',
    'UnmatchedParenthesesError',
    'InvalidBracketError',
    'EmptyFormulaError',
    'InvalidFormulaError',
    'UnsupportedAtomSpeciesError',
    'InvalidCountError',
]
