# Copyright 2026 Hiroshi Nakayama and RIKEN
# SPDX-License-Identifier: BSD-3-Clause

"""Peptide modification and b/y product-ion calculations.

Residue formulas are combined with a small public API to construct
b/y-type ion formulas and calculate their m/z values. The same
prefix/suffix approach can be adapted to other product-ion series.

Phosphorylation and H3PO4 loss are represented as reusable formula
changes. For illustration, the H3PO4-loss atom-species counts of pSer
are treated as a dehydroalanine residue in the conventional residue-list
fragmentation model used below.

This model does not represent interresidue connectivity changes
or reaction mechanisms. Gas-phase H3PO4 loss from phosphoserine may
involve a cyclic oxazoline formation with the neighboring N-terminal
carbonyl; such products require a connectivity-aware fragmentation model.
"""

from molformula import Formula, IonFormula, ion_mz, make_formula_change, protonated, sum_formulas


def b_y_ion_series(
        seq: list[Formula],
        z: int
        ) -> tuple[list[IonFormula], list[IonFormula]]:
    water = Formula('H2O')
    b_ions: list[IonFormula] = []
    y_ions: list[IonFormula] = []
    n = len(seq)
    for i in range(n):
        b_neutral = sum_formulas(*seq[: i+1])
        b_ion = protonated(b_neutral, z)
        b_ions.append(b_ion)

        y_neutral = sum_formulas(water, *seq[n-1-i: ])
        y_ion = protonated(y_neutral, z)
        y_ions.append(y_ion)

    return b_ions, y_ions


def mzs_for_products(
        products: list[IonFormula]) -> list[float]:
    return [ion_mz(p) for p in products]


def print_product_mzs(
        product_mzs: list[float],
        product_type: str
        ) -> None:
    for i, product_mz in enumerate(product_mzs):
        print(f'{product_type}{i+1}: {product_mz}')


if __name__ == '__main__':
    # Definition of amino acid residues
    ala = Formula('HNCH(CH3)CO')
    gly = Formula('HNCH2CO')
    ser = Formula('HNCH(CH2OH)CO')
    lys = Formula('HNCH(CH2CH2CH2CH2NH2)CO')

    phosphorylated = make_formula_change(
        gains=(Formula('HPO3'),)
    )
    h3po4_lost = make_formula_change(
        losses=(Formula('H3PO4'),)
    )
    pser = phosphorylated(ser)
    h3po4_lost_pser = h3po4_lost(pser)

    # Construct amino acid sequences
    unmodified_seq = [ala, gly, ser, lys]
    phosphorylated_seq = [ala, gly, pser, lys]
    # Treat the H3PO4-loss atom-species counts as a Dha residue
    # for the conventional b/y fragmentation model used below.
    h3po4_lost_seq = [ala, gly, h3po4_lost_pser, lys]

    # Product ion m/z calculations and outputs
    labeled_sequences = [
        ('Unmodified AGSK', unmodified_seq),
        ('AGpSK (Ser3 phosphorylated)', phosphorylated_seq),
        ('AGpSK after H3PO4 loss (hypothetical Dha residue)', h3po4_lost_seq)
    ]

    for label, sequence in labeled_sequences:
        b_ions, y_ions = b_y_ion_series(sequence, 1)
        b_ion_mzs = mzs_for_products(b_ions)
        y_ion_mzs = mzs_for_products(y_ions)

        print(f'Product ions for {label}')
        print_product_mzs(b_ion_mzs, 'b')
        print_product_mzs(y_ion_mzs, 'y')
        print()
