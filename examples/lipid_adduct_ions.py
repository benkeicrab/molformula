# Copyright 2026 Hiroshi Nakayama and RIKEN
# SPDX-License-Identifier: BSD-3-Clause

"""Operational construction of lipid adduct ions using composable formula changes.

The examples construct [M+H]+, [M+Na]+, and [M+HCOO]- by separating
atom-species-count changes from ionization:

                 formula change     ionization
------------------------------------------------------
    [M+H]+    :  -                  protonation
    [M+Na]+   :  +Na -H             protonation
    [M+HCOO]- :  +HCOOH             deprotonation

These operations construct ion formulas by composing atom-species-count
changes with ionization. They do not represent adduct-formation mechanisms
or predict whether an adduct is formed or observed experimentally.
"""

from typing import Callable

from molformula import Formula, IonFormula, apply_formula_change, deprotonated, ion_mz, make_formula_change, protonated, to_str


# Construct representative lipid formulas from building blocks via dehydration
h2o = Formula('H2O')
glycerol = Formula('C3H8O3')
phosphocholine = Formula('C5H14NO4P')
sphingosine = Formula('C18H37NO2')
palmitic_acid = Formula('C16H32O2')
oleic_acid = Formula('C18H34O2')

pc_34_1 = apply_formula_change(
    glycerol,
    gains=(phosphocholine, palmitic_acid, oleic_acid),
    losses=(h2o, h2o, h2o)
)

sm_34_1 = apply_formula_change(
    sphingosine,
    gains=(phosphocholine, palmitic_acid),
    losses=(h2o, h2o)
)

cer_34_1 = apply_formula_change(
    sphingosine,
    gains=(palmitic_acid,),
    losses=(h2o,)
)

lipids = [pc_34_1, sm_34_1, cer_34_1]

# Reusable atom-species-count changes for adduct-ion construction
identity = make_formula_change()  # identity mapping
sodium_substituted = make_formula_change(
    gains=(Formula('Na'),),
    losses=(Formula('H'),),
)
formic_acid_added = make_formula_change(
    gains=(Formula('HCOOH'),)
)

# Compose adduct-ion generation pipelines and apply across the lipid list
z = 1
protonated_mzs = [ion_mz(protonated(identity(l), z))
                  for l in lipids]
print(f'[M+H]+:    {protonated_mzs}')
sodium_adduct_mzs = [ion_mz(protonated(sodium_substituted(l), z))
                     for l in lipids]
print(f'[M+Na]+:   {sodium_adduct_mzs}')
formate_adduct_mzs = [ion_mz(deprotonated(formic_acid_added(l), z))
                      for l in lipids]
print(f'[M+HCOO]-: {formate_adduct_mzs}')


# Optional: The same calculations can be expressed as reusable transformations
# by composing formula changes with ionization functions.
def compose_ionization(
        formula_change: Callable[[Formula], Formula],
        ionization: Callable[[Formula, int], IonFormula],
        z: int,
        ) -> Callable[[Formula], IonFormula]:
    def _composed(f: Formula) -> IonFormula:
        return ionization(formula_change(f), z)

    return _composed


z = 1
ion_transformations = [
    ('[M+H]+', compose_ionization(identity, protonated, z)),
    ('[M+Na]+', compose_ionization(sodium_substituted, protonated, z)),
    ('[M+HCOO]-', compose_ionization(formic_acid_added, deprotonated, z))
]
for lipid in lipids:
    for label, transform in ion_transformations:
        ion = transform(lipid)
        print(f'{to_str(lipid)}, {label}: '
              f'{to_str(ion.formula)} / {ion_mz(ion)}')