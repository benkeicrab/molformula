"""Public contracts for factory-only IonFormula creation.

This module uses only names exported from molformula. It does not call
IonFormula._from_parts or access the private charge exception.
"""
from collections.abc import Callable
from dataclasses import FrozenInstanceError
from decimal import Decimal
from typing import Any

import pytest

from molformula import Formula, EmptyIonFormulaError, InsufficientHydrogenError, IonFormula, InvalidProtonCountError, atom_species_counts, deelectronated, deprotonated, ion_mz, protonated


def _protonated(formula: Formula) -> IonFormula:
    return protonated(formula, 2)


def _deprotonated(formula: Formula) -> IonFormula:
    return deprotonated(formula, 2)


CREATORS = [
    pytest.param(_protonated, "C6H14O6", 2, id="protonated"),
    pytest.param(_deprotonated, "C6H10O6", -2, id="deprotonated"),
    pytest.param(deelectronated, "C6H12O6", 1, id="deelectronated"),
]


@pytest.mark.parametrize("style", ["empty", "positional", "keyword"])
def test_direct_construction_is_rejected(style: str) -> None:
    with pytest.raises(TypeError):
        if style == "empty":
            IonFormula()
        elif style == "positional":
            IonFormula(Formula("H"), 1)
        else:
            IonFormula(formula=Formula("H1"), z=1)


@pytest.mark.parametrize("create, expected, z", CREATORS)
def test_factory_preserves_input_and_sets_ion_values(
    create: Callable[[Formula], IonFormula], expected: str, z: int,
) -> None:
    source = Formula("C6H12O6")
    before = atom_species_counts(source)
    result = create(source)
    assert isinstance(result, IonFormula)
    assert result.formula == Formula(expected)
    assert result.z == z
    assert atom_species_counts(source) == before


@pytest.mark.parametrize("create, expected, z", CREATORS)
def test_factory_results_have_value_equality_and_hash(
    create: Callable[[Formula], IonFormula], expected: str, z: int,
) -> None:
    first = create(Formula("C6H12O6"))
    second = create(Formula("C6H12O6"))
    assert first == second
    assert hash(first) == hash(second)


@pytest.mark.parametrize("create, expected, z", CREATORS)
@pytest.mark.parametrize("attribute", ["formula", "z"])
def test_factory_results_are_immutable(
    create: Callable[[Formula], IonFormula], expected: str,
    z: int, attribute: str,
) -> None:
    ion = create(Formula("C6H12O6"))
    replacements: dict[str, Any] = {
        "formula": Formula("C"), "z": 0,
    }
    with pytest.raises(FrozenInstanceError):
        setattr(ion, attribute, replacements[attribute])


@pytest.mark.parametrize("text", ["C6H12O6", "H", "O", "[2H]2O", "C[13C]N"])
def test_deelectronated_does_not_change_atom_counts(text: str) -> None:
    source = Formula(text)
    result = deelectronated(source)
    assert result.formula == source
    assert atom_species_counts(result.formula) == atom_species_counts(source)
    assert result.z == 1


@pytest.mark.parametrize("create, expected, z", CREATORS)
@pytest.mark.parametrize("bad", [None, "C6H12O6", 1])
def test_factories_reject_non_formula_input(
    create: Callable[[Formula], IonFormula], expected: str,
    z: int, bad: Any,
) -> None:
    with pytest.raises(TypeError):
        create(bad)


def test_deelectronated_rejects_an_ion_as_input() -> None:
    ion = protonated(Formula("H2O"), 1)
    with pytest.raises(TypeError):
        deelectronated(ion)  # type: ignore[arg-type]


def test_hydrogen_absence_is_allowed_for_protonation() -> None:
    result = protonated(Formula("O"), 1)
    assert result.formula == Formula("HO")
    assert result.z == 1


def test_removing_all_hydrogen_is_allowed_if_other_atoms_remain() -> None:
    result = deprotonated(Formula("HO"), 1)
    assert result.formula == Formula("O")
    assert result.z == -1


def test_hydrogen_shortage_has_public_error_and_attributes() -> None:
    with pytest.raises(InsufficientHydrogenError) as caught:
        deprotonated(Formula("H"), 2)
    assert caught.value.requested == 2
    assert caught.value.available == 1


def test_empty_ion_has_distinct_public_error() -> None:
    with pytest.raises(EmptyIonFormulaError):
        deprotonated(Formula("H"), 1)


@pytest.mark.parametrize("create", [protonated, deprotonated])
@pytest.mark.parametrize("bad", [True, False, 1.5])
def test_proton_counts_reject_bool_and_float(create: Callable[..., IonFormula], bad: Any) -> None:
    with pytest.raises(TypeError):
        create(Formula("C6H12O6"), bad)


@pytest.mark.parametrize("create", [protonated, deprotonated])
@pytest.mark.parametrize("bad", [0, -1])
def test_proton_counts_reject_nonpositive_values(create: Callable[..., IonFormula], bad: int) -> None:
    with pytest.raises(InvalidProtonCountError) as caught:
        create(Formula("C6H12O6"), bad)
    assert caught.value.n_proton == bad


@pytest.mark.parametrize("create, expected, z", CREATORS)
def test_ion_mz_matches_final_formula_and_charge(
    create: Callable[[Formula], IonFormula], expected: str, z: int,
) -> None:
    # Independent arithmetic with the decimal constants in the supplied source;
    # this checks electron correction, not external accuracy of those constants.
    masses = {"C": Decimal("12"), "H": Decimal("1.00782503223"),
              "O": Decimal("15.99491461957")}
    counts = {2: {"C": 6, "H": 14, "O": 6},
              -2: {"C": 6, "H": 10, "O": 6},
              1: {"C": 6, "H": 12, "O": 6}}[z]
    formula_mass = sum((masses[symbol] * count for symbol, count in counts.items()), Decimal(0))
    electron = Decimal("0.0005485799090441")
    expected_mz = float((formula_mass - z * electron) / abs(z))
    result = create(Formula("C6H12O6"))
    assert ion_mz(result) == pytest.approx(expected_mz, rel=0, abs=1e-10)
