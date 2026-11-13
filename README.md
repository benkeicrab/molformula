# molformula

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24096461.svg)](https://doi.org/10.5281/zenodo.24096461)
[![ChemRxiv](https://img.shields.io/badge/ChemRxiv-10.26434%2Fchemrxiv.15012007%2Fv1-blue)](https://doi.org/10.26434/chemrxiv.15012007/v1)

A simple, dependency-free Python library for safe chemical-formula manipulation and mass-spectrometry-oriented calculations.

This package is a reference implementation of the architecture described in our preprint:
> Hiroshi Nakayama and Masami Koike, **"Composable Chemical Formula Transformations with Invariant-Preserving Representations for Mass Spectrometry,"** *ChemRxiv* (2026). https://doi.org/10.26434/chemrxiv.15012007/v1

## Overview

This library handles molecular formulas and simple condensed formulas (See Supported formula syntax), arithmetic processing of these formulas, ionization, and exact mass and m/z calculations with enforced formula invariants.
The API is designed to support a functional style, using immutable formula values and composable, non-mutating operations.
It is intended as a foundation for mass spectrometry data processing.

Interactive demo: [![JupyterLite](https://jupyterlite.rtfd.io/en/latest/_static/badge.svg)](https://benkeicrab.github.io/molformula/lab/index.html?path=demo.ipynb)


## Installation

### Python version

Supports Python >=3.10

### Dependency

No third-party library is required for using the library. `pytest` for developers / contributors.

### Installation / uninstallation

```Bash
# installation
pip install molformula

# uninstallation
pip uninstall molformula
```

## Usage

### Import

```python
>>> import molformula as mf
```
or more explicitly (and thus more preferable)

```python
>>> from molformula import Formula, IonFormula, apply_formula_change, deprotonated, exact_mass, ion_mz, make_formula_change, protonated, sum_formulas, to_str

```

### Creating and showing canonical formulas

```python
>>> # molecular formula and simple condensed formulas are supported.
>>> glucose = Formula("CH2OH(CHOH)4CHO")
>>> # Output is canonicalized to Hill-system order [1].
>>> to_str(glucose)
'C6H12O6'
>>> # A symbol for explicit isotope must be written with a bracket in a formula.
>>> labeled_glucose = Formula("[13C]H2OH(CHOH)4CHO")
>>> # Labeled positions in a condensed formula are not preserved.
>>> # this package's additional rule: the ordinary symbol first, then explicit isotopes in increasing mass-number order.
>>> to_str(labeled_glucose)
'C5[13C]H12O6'
```

### Formula arithmetic

```python
>>> H2O = Formula('H2O')
>>> apply_formula_change(glucose, losses=(H2O,))
Formula('C6H10O5')
>>> C, C13 = Formula('C'), Formula('[13C]')
>>> apply_formula_change(glucose, gains=(C13,), losses=(C,))
Formula('C5[13C]H12O6')
```

### Reusable formula changes

Repeating application of a chemical change is supported.

```python
>>> # Define a specific formula change such as dehydration.
>>> dehydrated = make_formula_change(losses=(H2O,))
>>> dehydrated(glucose)
Formula('C6H10O5')
>>> C13_labeled = make_formula_change(gains=(C13,), losses=(C,))
>>> C13_labeled(glucose)
Formula('C5[13C]H12O6')
>>> C13_labeled(C13_labeled(glucose))
Formula('C4[13C]2H12O6')
```

### Ionization

```python
>>> protonated_glucose = protonated(glucose, 1)
>>> protonated_glucose
IonFormula(formula=Formula('C6H13O6'), z=1)
>>> deprotonated_glucose = deprotonated(glucose, 1)
>>> deprotonated_glucose
IonFormula(formula=Formula('C6H11O6'), z=-1)
```

### Exact mass and m/z

```python
>>> exact_mass(glucose)
180.06338810418
>>> ion_mz(deprotonated_glucose)
179.05611165185906
```

Exact mass and m/z are calculated using atomic mass values from NIST Physical Measurement Laboratory [2] and the electron mass from the 2022 CODATA recommendations [3].


### Application in peptide MS/MS

The following example constructs b- and y-ion series from amino-acid residue formulas for peptide AGSK.

```python
>>> def b_y_ion_series(seq, z):
...     water = Formula('H2O')
...     b_ions, y_ions = [], []
...     n = len(seq)
...     for i in range(n):
...         b_neutral = sum_formulas(*seq[: i+1])
...         b_ion = protonated(b_neutral, z)
...         b_ions.append(b_ion)
...
...         y_neutral = sum_formulas(water, *seq[n-1-i: ])
...         y_ion = protonated(y_neutral, z)
...         y_ions.append(y_ion)
...
...     return b_ions, y_ions
...
>>> gly = Formula('HNCH2CO')
>>> ser = Formula('HNCH(CH2OH)CO')
>>> ala = Formula('HNCH(CH3)CO')
>>> lys = Formula('HNCH(CH2CH2CH2CH2NH2)CO')
>>> pep_seq = [ala, gly, ser, lys]

>>> b_ions, y_ions = b_y_ion_series(pep_seq, 1)
>>> [ion_mz(ion) for ion in b_ions]
[72.04439023747095, 129.06585395816097, 216.09788236288097, 344.19284537807096]
>>> [ion_mz(ion) for ion in y_ions]
[147.11280415154096, 234.14483255626095, 291.16629627695096, 362.20341006210094]
```

The example includes full-length prefix and suffix compositions for completeness; these entries do not correspond to distinct backbone cleavage sites.

For phosphorylation of Ser and H3PO4-loss of phospho-Ser,
see examples/peptide_msms.py.


## Supported formula syntax

### Supported notation

- 'C6H12O6'        # Molecular formula
- 'C5[13C]H12O6'   # Molecular formula with explicit isotope
- 'CH4O'           # Count 1 can be omitted.
- 'CH3OH'          # Simple condensed formula
- '(CH3)2CHOH'     # Condensed formula with parentheses and multipliers
- '[ 13 C] 6H12O6' # Whitespace is allowed between tokens but not within lexical tokens such as elemental symbols or integer counts.

### Unsupported notation

- '2H2O'           # leading stoichiometric coefficient
- 'CuSO4.H2O'      # component separator / hydrate notation
- '[CH3COO]-'      # Charged formula
- 'K3[Fe(CN)6]'    # Complex
- 'CH3-OH'         # Bond representation with '-' and/or '='

- Formulas with unsupported atom species (See Limitations below)


## API overview

The library has immutable types: Formula and IonFormula, as well as functions for manipulating these types. 

### Types

#### `Formula`:

- Representation of non-negative atom-species counts.

- Must contain at least one atom.

- Does not encode charge, bonding, valence, or molecular structure.


#### `IonFormula`:

- Combines a Formula with a non-zero charge (z).

- Cannot be constructed directly; use ionization functions below.

- Ionization functions treat their Formula argument as a neutral
precursor.


### Operations

#### Parsing and formatting

- `Formula()`: Construct formula object.

- `to_str()`: Output a canonical representation of formula in Hill-system order [1]. Structural information in the original notation is not preserved.

- `atom_species_counts()`: Return atom-species counts as a new dictionary.

#### Arithmetic operations on Formulas

The library offers a high-level API of arithmetic operations on Formulas. 

Intermediates can have negative atom counts; however, the resulting Formula must have non-negative counts and contain at least one atom.

- `apply_formula_change()`:
    Apply atom-count gains and losses to a target Formula. 

- `make_formula_change()`:
    The atom-count gains and losses are baked into the return function. 
    The function applies gains/losses to a target Formula.

- `sum_formulas()`:
    Sum up multiple Formulas.


#### Ionization

Ionization functions manipulate atom counts and charges. Currently, three types of ionization are supported.

- `protonated()`:
    Adds n hydrogen atoms and sets charge to +n.

- `deprotonated()`:
    Removes n hydrogen atoms and sets charge to -n.

- `deelectronated()`:
    Leaves atom counts unchanged and sets charge to +1.


#### Mass calculation

`exact_mass()` calculates the exact mass of a `Formula`, and `ion_mz()`
calculates the m/z of an `IonFormula`.
The mass and m/z values are calculated with the atomic/isotopic masses [2] and electron mass [3].
No explicit decimal rounding is applied to the calculated mass or m/z values.

- `exact_mass()`:
    Calculates the exact mass of a Formula. No electron-mass correction is applied.

- `ion_mz()`:
    Calculates the ion m/z of an IonFormula.
    The m/z value of ion: IonFormula is calculated as

    (exact_mass(ion.formula) - ion.z * electron_mass) / abs(ion.z)

    so positive charge removes electron mass and negative charge adds it.


#### Exceptions

The class hierarchy is as shown below:

```bash
TypeError  # Error in Python types
ValueError
└── ChemicalFormulaError  # Base class for domain-value errors.
    ├── FormulaChangeError  # Change producing invalid Formula.
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
```

- The `position` attribute in `FormulaParseError` is a zero-based offset into the original input string. It may equal `len(input)` when expected input is missing at the end of the string.

- When raised by this package, the `deficits` attribute in `FormulaChangeError` is a dictionary mapping each atom species that would have a negative count to the positive number of missing atoms.

- The `zero_result` attribute in `FormulaChangeError` is True only when the formula change would remove all atoms exactly, without producing any negative counts.

- See the API documentation for details.


## Limitations

- Supported atom species: C, 13C, H, 2H, N, 15N, O, 18O, P, S, Se, F, Cl, 37Cl, Br, 81Br, I, Si, Li, Na, K, Mg, Ca, Fe;
The principal isotope of each element is represented by the ordinary elemental symbol. Other supported isotopes use mass-number-prefixed atom-species symbols, which are written in brackets in formula strings (e.g. 13C as [13C]).

- Complexes, hydrates, SMILES, and condensed-formula syntax other than the supported forms shown above are not supported.

- Formula validates atom-species counts, but does not validate chemical valence, bonding, structural plausibility, or electrical neutrality.

- Arithmetic operations on IonFormula are not allowed; ion reactions are expressed by arithmetic operations on Formula followed by ionization.

- Radical character and spin state are not represented.

- `to_str()` returns a canonical formula in Hill-system order, losing input order, grouping, etc.

- Adduct ions are not supported as a dedicated model. Existing formula-change and ionization operations can be composed to construct many adduct compositions, for example:

```Python
>>> glycine = Formula('C2H5NO2')
>>> ammonium_adduct = make_formula_change(gains=(Formula('NH3'),))
>>> ammonium_adducted = protonated(ammonium_adduct(glycine), 1)
>>> ammonium_adducted
IonFormula(formula=Formula('C2H9N2O2'), z=1)
```

See examples/lipid_adduct_ions.py for more details.

However, molformula does not determine whether such ions are chemically plausible or experimentally observable.

- `deprotonated()` strictly removes ordinary hydrogen (1H) atoms. It does not remove deuterium (2H) or tritium. To model deuterated loss, use `apply_formula_change()` with an explicit deuterium Formula('[2H]') as a loss.
```Python
>>> deuterium = Formula('[2H]')
>>> deuterated_glycine = Formula('C2H8[2H]N2O2')
>>> to_str(apply_formula_change(deuterated_glycine, losses=(deuterium,)))
'C2H8N2O2'
```

- Exact masses are calculated from the configured principal or explicitly specified isotopic masses. Isotopic abundance distributions, fine structure, and peak centroiding are not supported.


## License

This project is licensed under the 3-clause BSD License. See LICENSE for details.

## Citation

If you use this software, please cite the preprint above.


## Contributing

Bug reports and suggestions are welcome through GitHub Issues.
Contribution guidelines will be added as the project develops.


## References
1. Hill, E. A. On a system of indexing chemical literature; adopted by the classification division of the U. S. patent office. J. Am. Chem. Soc. 1900, 22 (8): 478–494.  https://doi.org/10.1021/ja02046a005

2. Atomic mass from NIST Physical Measurement Laboratory  https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl

3. 2022 CODATA recommended value  https://physics.nist.gov/cgi-bin/cuu/Value?eqmeu|search_for=electron+mass
