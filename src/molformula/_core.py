# Copyright 2026 Hiroshi Nakayama and RIKEN
# SPDX-License-Identifier: BSD-3-Clause

"""_core -- Chemical formula types and functions

This library aims to safely handle chemical formulas originally
developed for analyzing mass spectrometry data. The library has
Formula and IonFormula types, as well as functions for manipulating
these types.

Formula is an immutable value representing non-negative atom-species
counts and must contain at least one atom. It does not encode charge,
bonding, valence, or molecular structure.
The library offers high-level API of arithmetic operations of formulas
(apply_formula_change, make_formula_change, sum_formulas).

IonFormula combines a Formula with a non-zero charge.
Ionization functions treat their Formula argument as a neutral
precursor where documented.

Formula and IonFormula are transformed to exact mass
by exact_mass and m/z by ion_mz, respectively.

Types:
    Formula
        Public validated formula

    _FormulaRep
        Private fixed-width representation
        Signed counts are allowed for intermediate calculations

    IonFormula
        Ion result with Formula and non-zero charge

Public functions:
1. I/O
    Formula
        Parses a formula string and returns an instance of Formula

    to_str
        Returns a canonical representation of formula string

2. Arithmetic operations
    apply_formula_change
        Apply atom-count gains and losses to a target Formula.

    make_formula_change
        The atom-count gains and losses are baked into the return function
        which applies gains/losses to a target Formula.

    sum_formulas
        Sum up multiple Formulas.

3. Ionization
    protonated
        Increases hydrogen counts and charges

    deprotonated
        Decreases hydrogen counts and charges

    deelectronated
        Atom counts are not changed, singly positive charged

4. Mass calculation
    exact_mass
        Calculates exact mass of Formula

    ion_mz
        Calculates an ion m/z from an instance of IonFormula


"""

from dataclasses import dataclass
from enum import Enum
from math import fabs
from types import MappingProxyType
from typing import Callable, Final, Iterable, TypeAlias

import molformula._error as _error
import molformula._parser as _parser
from molformula._hill import _AtomKey, _format_hill


class _AtomSpecies(Enum):
    """Supported atom species.

    Each value is:
        (elemental_symbol, explicit_isotope_mass_number, exact_mass)

    None for the mass number denotes the ordinary elemental symbol,
    whose mass is that of the principal isotope.

    Atomic mass from NIST Physical Measurement Laboratory,
    retrieved 2026-10-05,
    https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=html&isotype=some
    """
    C = ('C', None, 12.00000000)
    C13 = ('C', 13, 13.003_354_835_07)
    H = ('H', None, 1.007_825_032_23)
    H2 = ('H', 2, 2.014_101_778_12)
    N = ('N', None, 14.003_074_004_43)
    N15 = ('N', 15, 15.000_108_898_88)
    O = ('O', None, 15.994_914_619_57)
    O18 = ('O', 18, 17.999_159_612_86)
    P = ('P', None, 30.973_761_998_42)
    S = ('S', None, 31.972_071_174_4)
    Se = ('Se', None, 79.916_521_8)
    F = ('F', None, 18.998_403_162_73)
    Cl = ('Cl', None, 34.968_852_682)
    Cl37 = ('Cl', 37, 36.965_902_602)
    Br = ('Br', None, 78.918_337_6)
    Br81 = ('Br', 81, 80.916_289_7)
    I = ('I', None, 126.904_471_9)
    Si = ('Si', None, 27.976_926_534_65)
    Li = ('Li', None, 7.016_003_436_6)
    Na = ('Na', None, 22.989_769_282_0)
    K = ('K', None, 38.963_706_4864)
    Mg = ('Mg', None, 23.985_041_697)
    Ca = ('Ca', None, 39.962_590_863)
    Fe = ('Fe', None, 55.934_936_33)


def _symbol(e: _AtomSpecies) -> str:
    elemental_symbol, mass_number, _ = e.value
    return (
        elemental_symbol
        if mass_number is None
        else f"{mass_number}{elemental_symbol}"
    )


_ATOM_SPECIES_SYMBOLS = tuple(_symbol(e) for e in _AtomSpecies)
_ATOM_SPECIES_SYMBOLS_SET = frozenset(_ATOM_SPECIES_SYMBOLS)
_ATOM_SPECIES_MASSES = tuple(e.value[2] for e in _AtomSpecies)
_INDEX_BY_SYMBOL: Final = MappingProxyType({_symbol(e): i
                                            for i, e in enumerate(_AtomSpecies)})
_ATOM_KEY: dict[str, _AtomKey] = {_symbol(e): (e.value[0], e.value[1])
                                  for e in _AtomSpecies}

# Electron mass in unified atomic mass units (u).
# 2022 CODATA recommended value, retrieved 2026-10-05.
# https://physics.nist.gov/cgi-bin/cuu/Value?eqmeu|search_for=electron+mass
_MASS_ELECTRON = 5.485_799_090_441e-4

FormulaStr: TypeAlias = str
Charge: TypeAlias = int


@dataclass(frozen=True, slots=True)
class _FormulaRep:
    """ Internal representation of chemical formula
    atom species counts can be zero and negative integer for intermediates.
    The available atom species are defined by enum _AtomSpecies.
    """
    _counts: tuple[int, ...]

    def __post_init__(self) -> None:
        if type(self._counts) is not tuple:
            raise TypeError("Counts must be a tuple.")

        if len(self._counts) != len(_AtomSpecies):
            raise ValueError(
                f"expected {len(_AtomSpecies)} atom-species counts, "
                f"got {len(self._counts)}"
            )

        if any(type(n) is not int for n in self._counts):
            raise TypeError("Atom-species counts must be integers.")


@dataclass(frozen=True, slots=True, init=False, repr=False, match_args=False,)
class Formula:
    """Immutable chemical formula represented by atom-species counts.

    Atom counts are non-negative and at least one atom must be present.
    Input ordering, grouping, and condensed-formula structure are not
    preserved. Charge, bonding, and valence are not represented.

    Equivalent formula strings therefore produce equal Formula values;
    for example, ``Formula("CH3OH") == Formula("CH4O")``.
    """
    _inner: _FormulaRep

    def __init__(self, formula_str: FormulaStr) -> None:
        """Construct a Formula from chemical formula string.

        Examples:
        >>> to_str(Formula('C5[13C]H12O6'))
        'C5[13C]H12O6'
        """
        rep = _rep_from(formula_str)
        _validate_rep(rep)
        object.__setattr__(self, "_inner", rep)

    @classmethod
    def _from_rep(cls, u: _FormulaRep) -> "Formula":
        """Internal factory for the calculation results
        of _add, _addall, _sub, and _mul
        """
        obj = object.__new__(cls)
        object.__setattr__(obj, "_inner", u)
        return obj

    def __str__(self) -> str:
        return to_str(self)

    def __repr__(self) -> str:
        return f"Formula({to_str(self)!r})"


@dataclass(frozen=True, slots=True, init=False, match_args=False,)
class IonFormula:
    """Immutable representation of an ion.

    Contains a Formula and a non-zero charge. Equality
    includes formula, and z. IonFormula cannot be constructed
    directly; use protonated(), deprotonated(), or deelectronated().

    Radical character and spin state are not represented.

    Examples:
    >>> protonated(Formula('C6H12O6'), 2)
    IonFormula(formula=Formula('C6H14O6'), z=2)
    >>> deprotonated(Formula('C6H12O6'), 2)
    IonFormula(formula=Formula('C6H10O6'), z=-2)
    >>> deelectronated(Formula('C6H12O6'))
    IonFormula(formula=Formula('C6H12O6'), z=1)
    """
    formula: Formula
    z: Charge

    def __init__(self, *args: object, **kwargs: object) -> None:
        raise TypeError(
            "IonFormula cannot be constructed directly; "
            "use an ionization function."
        )

    @classmethod
    def _from_parts(
        cls,
        *,
        formula: Formula,
        z: Charge,
    ) -> "IonFormula":
        """Construct an ion from its final formula and charge.

        Examples:
        >>> IonFormula._from_parts(formula=Formula('C6H12O6'), z=2)
        IonFormula(formula=Formula('C6H12O6'), z=2)
        >>> IonFormula._from_parts(formula=Formula('C6H12O6'), z=-2)
        IonFormula(formula=Formula('C6H12O6'), z=-2)
        >>> IonFormula._from_parts(formula=Formula('C6H12O6'), z=1.5)
        Traceback (most recent call last):
            ...
        TypeError: Ion charge must be an integer.
        >>> IonFormula._from_parts(formula=Formula('H1'), z=0)
        Traceback (most recent call last):
            ...
        molformula._error._InvalidChargeError: Ion charge must be non-zero, but got 0.
        """
        if not isinstance(formula, Formula):
            raise TypeError("Formula must be Formula.")

        _validate_charge(z)

        obj = object.__new__(cls)
        object.__setattr__(obj, "formula", formula)
        object.__setattr__(obj, "z", z)
        return obj


def _validate_charge(z: object) -> Charge:
    """

    Raises:
        TypeError: If an argument has an invalid Python type.
        _InvalidChargeError: If z is zero.

    Examples:
    >>> _validate_charge(1)
    1
    >>> _validate_charge(-1)
    -1
    >>> _validate_charge(1.0)
    Traceback (most recent call last):
        ...
    TypeError: Ion charge must be an integer.
    >>> _validate_charge(0)
    Traceback (most recent call last):
        ...
    molformula._error._InvalidChargeError: Ion charge must be non-zero, but got 0.
    """
    if type(z) is not int:
        raise TypeError("Ion charge must be an integer.")

    if z == 0:
        raise _error._InvalidChargeError(z)

    return z


def _validate_rep(u: _FormulaRep) -> None:
    if not isinstance(u, _FormulaRep):
        raise TypeError("_inner must be _FormulaRep.")
    if any(c < 0 for c in u._counts):
        deficits = {_ATOM_SPECIES_SYMBOLS[i]: -num
                    for i, num in enumerate(u._counts)
                    if num < 0}
        raise _error.FormulaChangeError(deficits=deficits)
    if not any(c > 0 for c in u._counts):
        raise _error.FormulaChangeError(zero_result=True)


def _validate(u: _FormulaRep) -> Formula:
    _validate_rep(u)
    return Formula._from_rep(u)


def _as_rep(vf: Formula) -> _FormulaRep:
    if not isinstance(vf, Formula):
        raise TypeError("A formula must be Formula.")
    return vf._inner


def _rep_from_counts(
    counts_by_symbol: dict[str, int],
) -> _FormulaRep:
    counts = [0] * len(_AtomSpecies)

    for symbol, count in counts_by_symbol.items():
        counts[_INDEX_BY_SYMBOL[symbol]] = count

    return _FormulaRep(tuple(counts))


def _rep_from(formula_str: FormulaStr) -> _FormulaRep:
    """Returns internal formula representation

    Examples:
    >>> _rep_from('C5[13C]H12O6')
    _FormulaRep(_counts=(5, 1, 12, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
    >>> _rep_from('C6H12O6')._counts
    (6, 0, 12, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    """
    counts = _parser.parse_formula_str(
        _ATOM_SPECIES_SYMBOLS_SET,
        formula_str,
    )
    return _rep_from_counts(counts)


def to_str(f: Formula) -> FormulaStr:
    """Return the canonical formula string in Hill-system order.

    Input ordering, repeated atom terms, grouping, and whitespace are
    not preserved.

    Examples:
    >>> to_str(Formula('CH3 OH'))
    'CH4O'
    >>> to_str('C6H12O6')
    Traceback (most recent call last):
        ...
    TypeError: A formula must be Formula.
    >>> to_str(Formula('C5[13C]H12O6'))
    'C5[13C]H12O6'
    """
    return _format_hill({_ATOM_KEY[ak]: count for ak, count in atom_species_counts(f).items()})


def atom_species_counts(f: Formula) -> dict[str, int]:
    """Return positive atom-species counts as a new dictionary.

    Isotopes are represented by keys such as `"13C"` and `"2H"`.
    """
    return {_ATOM_SPECIES_SYMBOLS[i]: num
            for i, num in enumerate(_as_rep(f)._counts) if num > 0}


def _add(rep1: _FormulaRep, rep2: _FormulaRep) -> _FormulaRep:
    """Returns an internal representation with sums of atom-species counts of rep1 and rep2.

    Examples:
    >>> hexose = _rep_from('C6H12O6')
    >>> to_str(_validate(_add(hexose, hexose)))
    'C12H24O12'
    """
    return _FormulaRep(tuple(map(lambda r1, r2: r1 + r2, rep1._counts, rep2._counts)))


def _addall(reps: list[_FormulaRep]) -> _FormulaRep:
    """Utility function for summing up given internal formula representations.

    Examples:
    >>> hexose = _rep_from('C6H12O6')
    >>> to_str(_validate(_addall([hexose, hexose, hexose])))
    'C18H36O18'
    """
    sum_counts = [0] * len(_AtomSpecies)
    for exp in reps:
        for i, _ in enumerate(sum_counts):
            sum_counts[i] += exp._counts[i]
    return _FormulaRep(tuple(sum_counts))


def _sub(rep1: _FormulaRep, rep2: _FormulaRep) -> _FormulaRep:
    """Returns a subtraction of rep2 from rep1.
    Negative atom-species counts are allowed.

    Returns:
        The count of each atom may be zero or negative.

    Examples:
    >>> hexose = _rep_from('C6H12O6')
    >>> water = _rep_from('H2O')
    >>> to_str(_validate(_sub(hexose, water)))
    'C6H10O5'
    >>> to_str(_validate(_sub(hexose, hexose)))
    Traceback (most recent call last):
        ...
    molformula._error.FormulaChangeError: The empty formula is invalid.
    """
    return _FormulaRep(tuple(map(lambda r1, r2: r1 - r2, rep1._counts, rep2._counts)))


def _mul(rep: _FormulaRep, m: int) -> _FormulaRep:
    """Returns a _FormulaRep with m fold counts of the original _FormulaRep.
    A non-positive integer is allowed.

    Examples:
    >>> to_str(_validate(_mul(_rep_from('C6H12O6'), 2)))
    'C12H24O12'
    >>> _mul(_rep_from('C6H12O6'), 0)
    _FormulaRep(_counts=(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
    >>> to_str(_validate(_mul(_rep_from('C6H12O6'), 0)))
    Traceback (most recent call last):
        ...
    molformula._error.FormulaChangeError: The empty formula is invalid.
    """
    return _FormulaRep(tuple([n * m for n in rep._counts]))


def _mass_from_rep(rep: _FormulaRep) -> float:
    return sum([m * n for m, n in zip(_ATOM_SPECIES_MASSES, rep._counts)])


def exact_mass(neutral_formula: Formula) -> float:
    """Return the exact mass represented by the formula, in u (Da).

    The mass is calculated as the sum of the configured isotopic masses.
    No electron-mass correction is applied.
    """
    return _mass_from_rep(_as_rep(neutral_formula))


def _validate_proton_count(n_proton: object) -> int:
    if type(n_proton) is not int:
        raise TypeError("Number of protons must be an integer.")

    if n_proton <= 0:
        raise _error.InvalidProtonCountError(n_proton)

    return n_proton


def protonated(neutral_formula: Formula, n_proton: int) -> IonFormula:
    """Construct a protonated cation from a Formula assumed to be neutral.

    Adds `n_proton` hydrogen atoms and sets the resulting charge to
    `+n_proton`.

    Raises:
        TypeError: If an argument has an invalid Python type.
        InvalidProtonCountError: If n_proton is not positive.

    Examples:
    >>> protonated(Formula('C6H12O6'), 2)
    IonFormula(formula=Formula('C6H14O6'), z=2)
    >>> protonated(Formula('C6H12O6'), 1.5)
    Traceback (most recent call last):
        ...
    TypeError: Number of protons must be an integer.
    >>> protonated(Formula('H'), -2)
    Traceback (most recent call last):
        ...
    molformula._error.InvalidProtonCountError: Number of protons must be positive: -2.
    >>> protonated(Formula('C6H12O6'), 0)
    Traceback (most recent call last):
        ...
    molformula._error.InvalidProtonCountError: Number of protons must be positive: 0.
    """
    _validate_proton_count(n_proton)
    formula = _validate(_add(_as_rep(neutral_formula), _mul(_rep_from('H'), n_proton)))
    return IonFormula._from_parts(formula=formula, z=n_proton)


def _h_count(vf: Formula) -> int:
    return _as_rep(vf)._counts[_INDEX_BY_SYMBOL['H']]


def deprotonated(neutral_formula: Formula, n_proton: int) -> IonFormula:
    """Constructs a deprotonated anion from a Formula assumed to be neutral.

    Subtracts `n_proton` hydrogen atoms and sets the resulting charge to
    `-n_proton`. `n_proton` must be positive.

    Raises:
        TypeError: If an argument has an invalid Python type.
        InvalidProtonCountError: If `n_proton` is not positive.
        InsufficientHydrogenError: If fewer than `n_proton` hydrogen atoms
            are present.
        EmptyIonFormulaError: If removing the hydrogen atoms would leave
            an empty formula.

    Examples:
    >>> deprotonated(Formula('C6H12O6'), 2)
    IonFormula(formula=Formula('C6H10O6'), z=-2)
    >>> deprotonated(Formula("HO"), 1)
    IonFormula(formula=Formula('O'), z=-1)
    >>> deprotonated(Formula('C6H12O6'), 1.5)
    Traceback (most recent call last):
        ...
    TypeError: Number of protons must be an integer.
    >>> deprotonated(Formula('C6H12O6'), -2)
    Traceback (most recent call last):
        ...
    molformula._error.InvalidProtonCountError: Number of protons must be positive: -2.
    >>> deprotonated(Formula('C6H12O6'), 0)
    Traceback (most recent call last):
        ...
    molformula._error.InvalidProtonCountError: Number of protons must be positive: 0.
    >>> deprotonated(Formula('H'), 2)
    Traceback (most recent call last):
        ...
    molformula._error.InsufficientHydrogenError: Cannot remove 2 protons: formula contains only 1 H atoms.
    >>> deprotonated(Formula('H'), 1)
    Traceback (most recent call last):
        ...
    molformula._error.EmptyIonFormulaError: Removing 1 proton(s) from 'H' leaves an empty formula.
    """
    _validate_proton_count(n_proton)
    neutral_rep = _as_rep(neutral_formula)
    available = _h_count(neutral_formula)

    if n_proton > available:
        raise _error.InsufficientHydrogenError(
            requested=n_proton,
            available=available,
        )

    ion_rep = _sub(neutral_rep, _mul(_rep_from('H'), n_proton))
    if not any(ion_rep._counts):
        raise _error.EmptyIonFormulaError(
            requested=n_proton,
            formula_str=to_str(neutral_formula),
        )

    return IonFormula._from_parts(formula=_validate(ion_rep), z=-n_proton)


def deelectronated(neutral_formula: Formula) -> IonFormula:
    """Constructs a cation by removing one electron from a neutral precursor formula.

    The input formula is assumed to represent a neutral species.
    Atom-species counts are unchanged, and the resulting charge is +1.

    This operation can represent molecular radical cations formed
    from closed-shell neutral precursors, e. g.
    by electron ionization (EI; M + e → M⁺· + 2e).
    Radical character and spin state are not tracked.

    Raises:
        TypeError: If neutral_formula is not Formula.

    Examples:
    >>> deelectronated(Formula('C6H12O6'))
    IonFormula(formula=Formula('C6H12O6'), z=1)
    >>> deelectronated(Formula("O"))
    IonFormula(formula=Formula('O'), z=1)
    """
    return IonFormula._from_parts(formula=neutral_formula, z=1)


def ion_mz(ion: IonFormula) -> float:
    """Return the m/z value of an ion.

    The value is calculated as

    (exact_mass(ion.formula) - ion.z * electron_mass) / abs(ion.z)

    so positive charge removes electron mass and negative charge adds it.

    Raises:
        TypeError: If `ion` is not an IonFormula.

    Examples:
    >>> ion_mz(protonated(Formula('C6H12O6'), 1))
    181.07066455650096
    >>> ion_mz(protonated(Formula('C6H12O6'), 2))
    91.03897050441095
    >>> ion_mz(deprotonated(Formula('C6H12O6'), 1))
    179.05611165185906
    >>> ion_mz(deprotonated(Formula('C6H12O6'), 2))
    89.02441759976905
    >>> ion_mz(Formula("H"))
    Traceback (most recent call last):
        ...
    TypeError: Ion must be IonFormula.
    """
    if not isinstance(ion, IonFormula):
        raise TypeError("Ion must be IonFormula.")
    return (exact_mass(ion.formula) - ion.z * _MASS_ELECTRON) / fabs(ion.z)


def _validated_formulas(
        values: Iterable[Formula],
        argument_name: str,
        ) -> tuple[Formula, ...]:
    """
    """
    result = tuple(values)

    if any(not isinstance(v, Formula) for v in result):
        raise TypeError(
            f"{argument_name} must contain only Formula."
        )

    return result


def _make_change_rep(
    gains: tuple[Formula, ...],
    losses: tuple[Formula, ...],
) -> _FormulaRep:
    """
    """
    terms: list[_FormulaRep] = [
        _as_rep(formula)
        for formula in gains
    ]

    terms.extend(
        _mul(_as_rep(formula), -1)
        for formula in losses
    )

    return _addall(terms)


def _apply_change_rep(
    formula: Formula,
    change: _FormulaRep,
) -> Formula:
    """
    """
    if not isinstance(formula, Formula):
        raise TypeError("A formula must be Formula.")

    return _validate(
        _add(
            _as_rep(formula),
            change,
        )
    )


def apply_formula_change(
    formula: Formula,
    /,
    *,
    gains: Iterable[Formula] = (),
    losses: Iterable[Formula] = (),
) -> Formula:
    """Apply atom-count gains and losses to a Formula.

    All gains and losses are combined before the final result is validated.
    Signed intermediate counts are therefore allowed internally, while the
    returned Formula must be non-empty and contain no negative counts.

    Raises:
        TypeError: If the inputs are not Formula objects as required.
        FormulaChangeError: If the final formula is empty or contains
            negative atom counts.

    Examples:
    >>> water = Formula('H2O')
    >>> glucose = Formula('C6H12O6')
    >>> apply_formula_change(glucose, losses=(water,))
    Formula('C6H10O5')
    >>> apply_formula_change(glucose, losses=(water, water))
    Formula('C6H8O4')
    """
    gain_formulas = _validated_formulas(gains, "gains")
    loss_formulas = _validated_formulas(losses, "losses")

    change = _make_change_rep(
        gain_formulas,
        loss_formulas,
    )

    return _apply_change_rep(formula, change)


def make_formula_change(
    *,
    gains: Iterable[Formula] = (),
    losses: Iterable[Formula] = (),
) -> Callable[[Formula], Formula]:
    """Return a reusable function that applies fixed gains and losses.

    The `gains` and `losses` iterables are consumed when this function is
    called; subsequent changes to the original iterables do not affect the
    returned function.

    Convenient for a specific chemical operation,
    e. g. dehydration or substitution of 12C with 13C.

    Examples:
    >>> water = Formula('H2O')
    >>> dehydration = make_formula_change(losses=(water,))
    >>> glucose = Formula('C6H12O6')
    >>> dehydration(glucose)
    Formula('C6H10O5')
    >>> dehydration(dehydration(glucose))
    Formula('C6H8O4')
    >>> didehydration = make_formula_change(losses=(water, water))
    >>> didehydration(glucose)
    Formula('C6H8O4')
    """
    gain_formulas = _validated_formulas(gains, "gains")
    loss_formulas = _validated_formulas(losses, "losses")

    change = _make_change_rep(
        gain_formulas,
        loss_formulas,
    )

    def apply(formula: Formula) -> Formula:
        return _apply_change_rep(formula, change)

    return apply


def sum_formulas(
    first: Formula,
    /,
    *rest: Formula,
) -> Formula:
    """Return the atom-count sum of one or more Formula objects.

    Examples:
    >>> sum_formulas(Formula("H2"), Formula("O"))
    Formula('H2O')
    """
    return apply_formula_change(
        first,
        gains=rest,
    )
