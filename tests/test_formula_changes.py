"""Contract tests for the three public formula-arithmetic functions.

Tests deliberately import only the package-level API. 
FormulaChangeError.deficits reports positive missing atom counts.
"""

from collections.abc import Callable, Iterable, Iterator
from typing import Any

import pytest

from molformula import Formula, FormulaChangeError, apply_formula_change, atom_species_counts, make_formula_change, sum_formulas


def _apply_using_factory(
    formula: Formula,
    *,
    gains: Iterable[Formula] = (),
    losses: Iterable[Formula] = (),
) -> Formula:
    return make_formula_change(gains=gains, losses=losses)(formula)


APPLIERS = [
    pytest.param(apply_formula_change, id="direct"),
    pytest.param(_apply_using_factory, id="generated"),
]


@pytest.mark.parametrize(
    "texts, expected",
    [
        (("C6H12O6",), "C6H12O6"),
        (("C6H12O6", "C6H12O6"), "C12H24O12"),
        (("C", "H4", "O1"), "CH4O"),
        (("C", "[13C]"), "C[13C]"),
    ],
    ids=["one-formula", "two-glucoses", "three-parts", "isotopes"],
)
def test_sum_formulas(texts: tuple[str, ...], expected: str) -> None:
    formulas = tuple(Formula(text) for text in texts)
    assert sum_formulas(formulas[0], *formulas[1:]) == Formula(expected)


def test_sum_requires_at_least_one_formula() -> None:
    with pytest.raises(TypeError):
        sum_formulas()  # type: ignore[call-arg]


@pytest.mark.parametrize("bad_position", [0, 1], ids=["first", "later"])
def test_sum_rejects_non_formula_arguments(bad_position: int) -> None:
    values: list[Any] = [Formula("C1"), Formula("H1")]
    values[bad_position] = "H1"
    with pytest.raises(TypeError):
        sum_formulas(values[0], *values[1:])


@pytest.mark.parametrize("apply", APPLIERS)
@pytest.mark.parametrize(
    "base_text, gain_texts, loss_texts, expected_text",
    [
        ("C6H12O6", (), ("H2O",), "C6H10O5"),
        ("C", ("H4O",), (), "CH4O"),
        ("C2H4", ("[13C]",), ("C",), "C[13C]H4"),
        ("H", ("H2",), ("H2",), "H"),
        ("H", ("O",), ("H",), "O"),
        ("C", (), (), "C"),
        ("C", ("H2O",), ("H2O",), "C"),
        ("C", ("H4O", "H2"), ("H2O",), "CH4"),
    ],
    ids=[
        "dehydration", "gain-only", "c13-substitution",
        "loss-first-would-be-negative", "loss-first-would-be-empty",
        "identity", "net-zero-change", "multiple-terms",
    ],
)
def test_change_matches_final_net_formula(
    apply: Callable[..., Formula],
    base_text: str,
    gain_texts: tuple[str, ...],
    loss_texts: tuple[str, ...],
    expected_text: str,
) -> None:
    base = Formula(base_text)
    gains = tuple(Formula(text) for text in gain_texts)
    losses = tuple(Formula(text) for text in loss_texts)
    inputs = (base, *gains, *losses)
    before = [atom_species_counts(formula) for formula in inputs]

    result = apply(base, gains=gains, losses=losses)

    assert isinstance(result, Formula)
    assert result == Formula(expected_text)
    assert [atom_species_counts(formula) for formula in inputs] == before


@pytest.mark.parametrize("apply", APPLIERS)
def test_negative_result_reports_formula_change_error(
    apply: Callable[..., Formula],
) -> None:
    with pytest.raises(FormulaChangeError) as caught:
        apply(Formula("H"), losses=(Formula("H2O2"),))

    assert caught.type is FormulaChangeError
    assert caught.value.zero_result is False
    assert set(caught.value.deficits) == {"H", "O"}
    # The deficits must be positive integers, so that zero-result must be False.
    assert caught.value.deficits['H'] == 1
    assert caught.value.deficits['O'] == 2


@pytest.mark.parametrize("apply", APPLIERS)
def test_empty_result_reports_formula_change_error(
    apply: Callable[..., Formula],
) -> None:
    with pytest.raises(FormulaChangeError) as caught:
        apply(Formula("H"), losses=(Formula("H"),))

    assert caught.type is FormulaChangeError
    assert caught.value.zero_result is True
    assert caught.value.deficits == {}


@pytest.mark.parametrize("apply", APPLIERS)
@pytest.mark.parametrize("bad", [None, "H", 1], ids=["none", "string", "integer"])
def test_change_rejects_non_formula_base(
    apply: Callable[..., Formula], bad: object,
) -> None:
    with pytest.raises(TypeError):
        apply(bad)


@pytest.mark.parametrize("factory", [False, True], ids=["direct", "factory"])
@pytest.mark.parametrize("argument", ["gains", "losses"])
@pytest.mark.parametrize(
    "bad_values", [None, ("H",), (Formula("H"), None)],
    ids=["not-iterable", "string-element", "mixed-elements"],
)
def test_invalid_change_terms_are_rejected_at_construction(
    factory: bool, argument: str, bad_values: Any,
) -> None:
    kwargs: dict[str, Any] = {argument: bad_values}
    with pytest.raises(TypeError):
        if factory:
            make_formula_change(**kwargs)
        else:
            apply_formula_change(Formula("C1"), **kwargs)


def test_factory_snapshots_terms_at_creation() -> None:
    gains = [Formula("H")]
    losses: list[Formula] = []
    change = make_formula_change(gains=gains, losses=losses)
    gains[:] = [Formula("O")]
    losses.append(Formula("C"))

    assert change(Formula("C")) == Formula("CH")
    assert change(Formula("C")) == Formula("CH")


def test_factory_consumes_generators_once_at_creation() -> None:
    events: list[str] = []

    def gain_terms() -> Iterator[Formula]:
        events.append("gain")
        yield Formula("O")

    def loss_terms() -> Iterator[Formula]:
        events.append("loss")
        yield Formula("H")

    change = make_formula_change(gains=gain_terms(), losses=loss_terms())
    assert sorted(events) == ["gain", "loss"]
    for _ in range(2):
        assert change(Formula("CH")) == Formula("CO")
    assert sorted(events) == ["gain", "loss"]


def test_failed_application_does_not_damage_reusable_change() -> None:
    dehydrate = make_formula_change(losses=(Formula("H2O"),))
    with pytest.raises(FormulaChangeError):
        dehydrate(Formula("H1"))
    assert dehydrate(Formula("C6H12O6")) == Formula("C6H10O5")
    assert dehydrate(Formula("C6H12O6")) == Formula("C6H10O5")


@pytest.mark.parametrize("factory", [False, True], ids=["direct", "factory"])
@pytest.mark.parametrize("argument", ["gains", "losses"])
def test_external_iterator_error_is_not_reclassified(
    factory: bool, argument: str,
) -> None:
    original_error = ValueError("failure in the caller's iterator")

    def terms() -> Iterator[Formula]:
        yield Formula("H")
        raise original_error

    kwargs = {argument: terms()}
    with pytest.raises(ValueError) as caught:
        if factory:
            make_formula_change(**kwargs)
        else:
            apply_formula_change(Formula("C"), **kwargs)

    assert caught.value is original_error
