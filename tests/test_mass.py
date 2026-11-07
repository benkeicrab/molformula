"""Numerical mass contracts tested exclusively through molformula's public API.

References (checked 2026-10-05; no network access at test execution):
  NIST Atomic Weights and Isotopic Compositions, Relative Atomic Mass column:
  https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=html&isotype=some
  NIST 2022 CODATA, electron mass in u:
  https://physics.nist.gov/cuu/Constants/Table/allascii.txt

Reference values are independently recorded here, not imported from _core.
Decimal computes the reference sums; the implementation returns float.
These tests verify the adopted atom-mass-sum/electron-correction convention,
not physical uncertainty, binding energies, or instrument mass accuracy.

Only _formula() knows the current molecular- and condensed-formula syntax.

Python >= 3.10; pytest only. No production-code changes are required.
"""

from decimal import Decimal, localcontext
from math import isfinite
from typing import Any

import pytest

from molformula import Formula, apply_formula_change, deelectronated, deprotonated, exact_mass, ion_mz, make_formula_change, protonated, sum_formulas


ATOM_MASS_U = {
    "C": Decimal("12"),
    "13C": Decimal("13.00335483507"),
    "H": Decimal("1.00782503223"),
    "2H": Decimal("2.01410177812"),
    "N": Decimal("14.00307400443"),
    "15N": Decimal("15.00010889888"),
    "O": Decimal("15.99491461957"),
    "18O": Decimal("17.99915961286"),
    "P": Decimal("30.97376199842"),
    "S": Decimal("31.9720711744"),
    "Se": Decimal('79.9165218'),
    "F": Decimal("18.99840316273"),
    "Cl": Decimal("34.968852682"),
    "37Cl": Decimal("36.965902602"),
    "Br": Decimal("78.9183376"),
    "81Br": Decimal("80.9162897"),
    "I": Decimal("126.9044719"),
    "Si": Decimal("27.97692653465"),
    "Li": Decimal("7.0160034366"),
    "Na": Decimal("22.9897692820"),
    "K": Decimal("38.9637064864"),
    "Mg": Decimal("23.985041697"),
    "Ca": Decimal("39.962590863"),
    "Fe": Decimal("55.93493633"),
}
ELECTRON_MASS_U = Decimal("0.0005485799090441")

# Numerical regression tolerances, NOT physical uncertainty estimates.
# At mass 1.3e5, REL_TOL permits about 1.3e-9 in the returned mass units.
REL_TOL = 1e-14
ABS_TOL = 1e-12
# Subtraction of two computed masses can lose relative precision.
SHIFT_ABS_TOL = 1e-10

GLUCOSE = {"C": 6, "H": 12, "O": 6}
P_CONTAINING = {"C": 38, "H": 49, "N": 15, "O": 29, "P": 4}
MIXED_ISOTOPES = {
    "C": 3, "13C": 2, "H": 5, "2H": 4, "N": 2,
    "15N": 1, "O": 3, "18O": 2, "P": 1, "S": 1,
}
LARGE_FORMULA = {symbol: count * 100 for symbol, count in P_CONTAINING.items()}


def _formula(counts: dict[str, int]) -> Formula:
    """Test-data adapter, not a new production API or alternate parser."""
    elems = []
    for symbol, count in counts.items():
        symbol = f'[{symbol}]' if symbol[0].isdigit() else symbol
        if count == 1:
            elems.append(f'{symbol}')
        elif count > 1:
            elems.append(f'{symbol}{count}')
        else:
            pass

    return Formula(''.join(elems))


def _reference_mass(counts: dict[str, int]) -> Decimal:
    """Use explicit test counts, never atom_species_counts(actual_result)."""
    with localcontext() as ctx:
        ctx.prec = 50
        return sum(
            (ATOM_MASS_U[symbol] * count for symbol, count in counts.items()),
            Decimal(0),
        )


def _reference_proton_transfer_mz(counts: dict[str, int], z: int) -> Decimal:
    # This is the package's chosen atomic-mass/electron correction convention.
    # Do not substitute a differently tabulated proton mass here.
    with localcontext() as ctx:
        ctx.prec = 50
        h_transfer = ATOM_MASS_U["H"] - ELECTRON_MASS_U
        return (_reference_mass(counts) + z * h_transfer) / abs(z)


def _assert_mass(actual: float, expected: Decimal | float) -> None:
    assert isinstance(actual, float)
    assert isfinite(actual)
    assert actual == pytest.approx(float(expected), rel=REL_TOL, abs=ABS_TOL)


@pytest.mark.parametrize("symbol", tuple(ATOM_MASS_U))
def test_single_atom_mass_matches_reference(symbol: str) -> None:
    _assert_mass(exact_mass(_formula({symbol: 1})), ATOM_MASS_U[symbol])


@pytest.mark.parametrize(
    "counts",
    [
        pytest.param({"H": 2, "O": 1}, id="water"),
        pytest.param(GLUCOSE, id="glucose"),
        pytest.param(P_CONTAINING, id="p-containing"),
        pytest.param({"H": 3, "P": 1, "O": 4}, id="phosphorus"),
        pytest.param({"C": 9, "H": 13, "N": 3, "O": 4, "S": 1}, id="sulfur"),
        pytest.param(MIXED_ISOTOPES, id="all-supported-isotopes"),
        pytest.param(LARGE_FORMULA, id="large-formula-not-a-sequence-model"),
    ],
)
def test_formula_mass_matches_independent_reference(counts: dict[str, int]) -> None:
    _assert_mass(exact_mass(_formula(counts)), _reference_mass(counts))


@pytest.mark.parametrize(
    "text, expected",
    [
        ("H2O", "18.01056468403"),
        ("C6H12O6", "180.06338810418"),
        ("C38H49N15O29P4", "1303.17710860693"),
    ],
)
def test_formula_exact_mass_reference_examples(text: str, expected: str) -> None:
    _assert_mass(exact_mass(Formula(text)), Decimal(expected))


@pytest.mark.parametrize(
    "light, heavy",
    [("C", "13C"), ("H", "2H"), ("N", "15N"), ("O", "18O")],
)
def test_isotope_replacement_mass_shift(light: str, heavy: str) -> None:
    base = _formula(P_CONTAINING)
    replace = make_formula_change(
        gains=(_formula({heavy: 1}),),
        losses=(_formula({light: 1}),),
    )
    shift = exact_mass(replace(base)) - exact_mass(base)
    expected = float(ATOM_MASS_U[heavy] - ATOM_MASS_U[light])
    assert shift == pytest.approx(expected, rel=0, abs=SHIFT_ABS_TOL)


@pytest.mark.parametrize(
    "loss_counts", [{"H": 2, "O": 1}, {"H": 3, "P": 1, "O": 4}],
    ids=["water", "phosphoric-acid"],
)
def test_neutral_loss_mass_difference(loss_counts: dict[str, int]) -> None:
    base = _formula(P_CONTAINING)
    product = apply_formula_change(base, losses=(_formula(loss_counts),))
    difference = exact_mass(base) - exact_mass(product)
    assert difference == pytest.approx(
        float(_reference_mass(loss_counts)), rel=0, abs=SHIFT_ABS_TOL,
    )


def test_mass_is_additive_for_formula_sum() -> None:
    # A complementary invariant, not a replacement for reference-value tests.
    first = _formula(GLUCOSE)
    second = _formula(P_CONTAINING)
    _assert_mass(
        exact_mass(sum_formulas(first, second)),
        exact_mass(first) + exact_mass(second),
    )


@pytest.mark.parametrize("z", [1, 2, 3, -1, -2, -3])
@pytest.mark.parametrize(
    "counts", [GLUCOSE, P_CONTAINING, MIXED_ISOTOPES],
    ids=["glucose", "p-containing", "mixed-isotopes"],
)
def test_proton_transfer_mz_matches_independent_reference(
    counts: dict[str, int], z: int,
) -> None:
    neutral = _formula(counts)
    create = protonated if z > 0 else deprotonated
    ion = create(neutral, abs(z))
    _assert_mass(ion_mz(ion), _reference_proton_transfer_mz(counts, z))


@pytest.mark.parametrize(
    "counts", [GLUCOSE, {"O": 1}, MIXED_ISOTOPES],
    ids=["glucose", "without-hydrogen", "mixed-isotopes"],
)
def test_deelectronated_mz_matches_independent_reference(counts: dict[str, int]) -> None:
    ion = deelectronated(_formula(counts))
    _assert_mass(ion_mz(ion), _reference_mass(counts) - ELECTRON_MASS_U)


@pytest.mark.parametrize("z", [1, 2, 3, -1, -2, -3])
def test_isotope_mz_spacing_is_divided_by_absolute_charge(z: int) -> None:
    create = protonated if z > 0 else deprotonated
    light = _formula(P_CONTAINING)
    heavy_counts = dict(P_CONTAINING)
    heavy_counts["C"] -= 1
    heavy_counts["13C"] = 1
    heavy = _formula(heavy_counts)
    light_ion = create(light, abs(z))
    heavy_ion = create(heavy, abs(z))
    difference = ion_mz(heavy_ion) - ion_mz(light_ion)
    expected = float((ATOM_MASS_U["13C"] - ATOM_MASS_U["C"]) / abs(z))
    assert difference == pytest.approx(expected, rel=0, abs=SHIFT_ABS_TOL)


@pytest.mark.parametrize("kind", ["protonated", "deprotonated", "deelectronated"])
def test_exact_mass_of_ion_formula_has_no_electron_correction(kind: str) -> None:
    neutral = _formula(GLUCOSE)
    final_counts = dict(GLUCOSE)
    if kind == "protonated":
        ion = protonated(neutral, 2)
        final_counts["H"] += 2
    elif kind == "deprotonated":
        ion = deprotonated(neutral, 2)
        final_counts["H"] -= 2
    else:
        ion = deelectronated(neutral)

    _assert_mass(exact_mass(ion.formula), _reference_mass(final_counts))
    # The charge correction belongs to ion_mz, not to Formula/exact_mass.
    assert exact_mass(ion.formula) - abs(ion.z) * ion_mz(ion) == pytest.approx(
        float(ion.z * ELECTRON_MASS_U), rel=0, abs=SHIFT_ABS_TOL,
    )


@pytest.mark.parametrize("bad", [None, "C1", 12, True])
def test_exact_mass_rejects_non_formula(bad: Any) -> None:
    with pytest.raises(TypeError):
        exact_mass(bad)


def test_exact_mass_rejects_ion_object() -> None:
    ion = deelectronated(_formula(GLUCOSE))
    with pytest.raises(TypeError):
        exact_mass(ion)  # type: ignore[arg-type]
