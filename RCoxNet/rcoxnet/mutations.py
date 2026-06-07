import pandas as pd

# variant types considered functional
FUNCTIONAL_VARIANTS = {
    'Missense_Mutation',
    'Nonsense_Mutation',
    'Frame_Shift_Ins',
    'Frame_Shift_Del',
    'Splice_Site',
    'In_Frame_Ins',
    'In_Frame_Del',
    'Translation_Start_Site',
}


def parse_mutations(mut_raw):
    """
    Keep only functional variants and build per-patient seed gene lists
    and per-gene mutation counts (used for RWR seed construction).
    """
    mut_func = mut_raw[
        mut_raw['Variant_Classification'].isin(FUNCTIONAL_VARIANTS)
    ].copy()

    seed_dict = (mut_func
                 .groupby('Tumor_Sample_Barcode')['Hugo_Symbol']
                 .apply(lambda x: sorted(x.unique().tolist()))
                 .to_dict())

    mut_count_dict = (mut_func
                      .groupby('Tumor_Sample_Barcode')['Hugo_Symbol']
                      .apply(lambda x: x.value_counts().to_dict())
                      .to_dict())

    print(f'Patients with functional mutations: {len(seed_dict)}')
    return seed_dict, mut_count_dict
