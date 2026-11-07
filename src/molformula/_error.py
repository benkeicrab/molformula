# Copyright 2026 Hiroshi Nakayama and RIKEN
# SPDX-License-Identifier: BSD-3-Clause

"""Exceptions representing public contracts.

FormulaChangeError:
    When raised by this package, the `deficits` attribute is
    a dictionary mapping each atom species that would have
    a negative count to the positive number of missing atoms.

    The `zero_result` attribute is True only when the formula change
    would remove all atoms exactly, without producing any negative counts.

FormulaParseError:
    Exceptions representing public formula-parsing contracts.

    The `position` attribute is a zero-based offset
    into the original input string. It may equal `len(input)`
    when expected input is missing at the end of the string.

The class hierarchy is as shown below:
ValueError
└── ChemicalFormulaError
    ├── FormulaChangeError
    ├── IonizationError
    │   ├── InvalidProtonCountError
    │   ├── InsufficientHydrogenError
    │   └── EmptyIonFormulaError
    └── FormulaParseError
        ├── EmptyParenthesisError
        ├── UnmatchedParenthesesError
        ├── InvalidBracketError
        ├── EmptyFormulaError
        ├── InvalidFormulaError
        ├── UnsupportedAtomSpeciesError
        └── InvalidCountError
"""


class ChemicalFormulaError(ValueError):
    """Base class for domain-value errors raised by this package.
    """


class FormulaChangeError(ChemicalFormulaError):
    """A formula change would produce an invalid Formula.
    """
    def __init__(
        self,
        *,
        deficits: dict[str, int] | None = None,
        zero_result: bool = False,
    ) -> None:
        self.deficits = deficits or {}
        self.zero_result = zero_result

        if zero_result:
            message = "The empty formula is invalid."
        else:
            message = (
                "Formula change requires unavailable atoms: "
                f"{self.deficits}."
            )

        super().__init__(message)


class IonizationError(ChemicalFormulaError):
    """Base class for ionization failures."""


class InvalidProtonCountError(IonizationError):
    """The requested proton number must be positive."""
    def __init__(
        self,
        n_proton: int,
    ) -> None:
        self.n_proton = n_proton

        super().__init__(
            "Number of protons must be positive: "
            f"{self.n_proton}."
        )


class _InvalidChargeError(IonizationError):
    """The requested charge is not allowed, such as z == 0."""
    def __init__(
        self,
        z: int,
    ) -> None:
        self.z = z

        super().__init__(
            f"Ion charge must be non-zero, but got {self.z}."
        )


class InsufficientHydrogenError(IonizationError):
    """The requested deprotonation exceeds available H atoms."""
    def __init__(
        self,
        requested: int,
        available: int,
    ) -> None:
        self.requested = requested
        self.available = available

        super().__init__(
            f"Cannot remove {requested} protons: "
            f"formula contains only {available} H atoms."
        )


class EmptyIonFormulaError(IonizationError):
    """Removing hydrogen(s) would leave an empty formula.
    """
    def __init__(
        self,
        requested: int,
        formula_str: str,
    ) -> None:
        self.requested = requested
        self.formula_str = formula_str

        super().__init__(
            f"Removing {requested} proton(s) "
            f"from '{formula_str}' leaves an empty formula."
        )


class FormulaParseError(ChemicalFormulaError):
    """Base class for formula parsing errors."""
    def __init__(self, message: str, position: int) -> None:
        self.position = position
        super().__init__(message)


class EmptyParenthesisError(FormulaParseError):
    """A parenthesized group is empty."""
    def __init__(self, position: int) -> None:
        super().__init__(
            "Empty parenthesis",
            position,
        )


class UnmatchedParenthesesError(FormulaParseError):
    """An opening or closing parenthesis is unmatched."""
    def __init__(self, position: int) -> None:
        super().__init__(
            "Unmatched parenthesis",
            position,
        )


class InvalidBracketError(FormulaParseError):
    """Bracketed isotope notation is malformed,
    including a missing or unexpected bracket, mass number,
    or elemental symbol.
    """
    def __init__(self, position: int) -> None:
        super().__init__(
            "Invalid bracket",
            position,
        )


class EmptyFormulaError(FormulaParseError):
    """The input formula is empty or contains only whitespace."""
    def __init__(self, position: int) -> None:
        super().__init__(
            "Empty formula",
            position,
        )


class InvalidFormulaError(FormulaParseError):
    """An atom species cannot be parsed."""
    def __init__(self, position: int) -> None:
        super().__init__(
            "Invalid formula",
            position,
        )


class UnsupportedAtomSpeciesError(FormulaParseError):
    """An atom species is not among the supported species.
    """
    def __init__(self, symbol: str, position: int) -> None:
        self.symbol = symbol
        super().__init__(
            f"Unsupported atom species {self.symbol}",
            position,
        )


class InvalidCountError(FormulaParseError):
    """An atom / group count is zero or negative.
    """
    def __init__(self, character: str, position: int) -> None:
        self.character = character
        super().__init__(
            f"Invalid count starts with {self.character}",
            position,
        )
