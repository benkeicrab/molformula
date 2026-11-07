"""Small contract tests for the package-level public API."""

import pytest

import molformula
from molformula import Formula, atom_species_counts, exact_mass, to_str


EXPECTED_PUBLIC_NAMES = {
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
}


def test_documented_public_names_are_exported() -> None:
    # A subset assertion allows intentional additions to the public API.
    assert EXPECTED_PUBLIC_NAMES == set(molformula.__all__)
    for name in EXPECTED_PUBLIC_NAMES:
        assert hasattr(molformula, name), f"missing package export: {name}"


def test_internal_aliases_are_not_package_exports() -> None:
    internal_names = {
        "FormulaStr", "Charge", "_FormulaRep", "_AtomSpecies", "_as_rep",
    }
    assert internal_names.isdisjoint(molformula.__all__)
    # Check only the root namespace, not access to private submodules.
    assert not hasattr(molformula, "FormulaStr")
    assert not hasattr(molformula, "Charge")


def test_atom_species_counts_returns_an_independent_dictionary() -> None:
    formula = Formula("C6H12O6")
    counts = atom_species_counts(formula)
    counts["C"] = 999
    counts["P"] = 10
    assert atom_species_counts(formula) == {"C": 6, "H": 12, "O": 6}


# both molecular formulas and condensed formulas
@pytest.mark.parametrize(
    "text",
    [
        "C6H12O6",
        "CH4O",
        "CH2[2H]2O",
        "C3H8O",
        "C5[13C]H12O6",
        "H3O4P"
    ],
)
def test_invariant_formula(text: str) -> None:
    assert to_str(Formula(text)) == text


# condensed formulas
@pytest.mark.parametrize(
    "text, expected",
    [
        ("CH3OH", "CH4O"),
        ("CH[2H]2OH", "CH2[2H]2O"),
        ("(CH3)2CHOH", "C3H8O"),
    ],
)
def test_condensed_formula_is_canonicalized(
        text: str, expected: str,
    ) -> None:
    assert to_str(Formula(text)) == expected


# condensed formulas
@pytest.mark.parametrize(
    "text",
    [
        "CH3OOH",
        "CH[2H]2OH",
        "(CH3)2CHOH",
        "H3PO4",
    ],
)
def test_parser_formatter(text: str) -> None:
    formula = Formula(text)
    assert formula == Formula(to_str(formula))


def test_formula_str_uses_canonical_format() -> None:
    assert str(Formula("CH3OH")) == "CH4O"
    assert str(Formula("H2NCHCH2SeHCOONa")) == "C3H6NNaO2Se"


def test_formula_repr_uses_canonical_format() -> None:
    assert repr(Formula("C[13C]H4")) == "Formula('C[13C]H4')"


def test_ordinary_explicit_isotope_pair() -> None:
    assert to_str(Formula('Ca[81Br]Br')) == 'Br[81Br]Ca'
    assert exact_mass(Formula('[81Br]Br')) == pytest.approx(78.9183376 + 80.9162897, rel=1e-14)
