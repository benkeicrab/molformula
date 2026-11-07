# Copyright 2026 Hiroshi Nakayama and RIKEN
# SPDX-License-Identifier: BSD-3-Clause

"""Atom-species count validation using acetaminophen.

The C2H2O (ketene) loss is consistent with reported gas-phase
MS/MS fragmentation of protonated acetaminophen.

Shiea et al., Rapid Commun. Mass Spectrom. 2020, 34, e8564. DOI: 10.1002/rcm.8564

CO2 loss is hypothetical and is used to illustrate validation.

This example evaluates formula changes, not reaction mechanisms
or the observability of product ions.
"""

from molformula import Formula, FormulaChangeError, make_formula_change, to_str


acetaminophen = Formula('C8H9NO2')

# Reusable formula changes
ketene_lost = make_formula_change(losses=(Formula('C2H2O'),))
co2_lost = make_formula_change(losses=(Formula('CO2'),))

# Reported neutral loss: valid atom-species counts
ketene_lost_acetaminophen = ketene_lost(acetaminophen)
assert to_str(ketene_lost_acetaminophen) == 'C6H7NO'
print(f'Ketene loss from {acetaminophen}: {ketene_lost_acetaminophen}')

# Hypothetical loss: valid atom-species counts, plausibility not assessed
co2_lost_acetaminophen = co2_lost(acetaminophen)
assert to_str(co2_lost_acetaminophen) == 'C7H9N'
print(f'Hypothetical CO2 loss from {acetaminophen}: {co2_lost_acetaminophen}')

# Invalid atom-species counts: one oxygen atom is unavailable
try:
    co2_lost(ketene_lost_acetaminophen)
except FormulaChangeError as err:
    assert err.deficits == {"O": 1}
    print(f'Rejected CO2 loss from {ketene_lost_acetaminophen}; deficits: {err.deficits}')
else:
    raise AssertionError("Expected FormulaChangeError")