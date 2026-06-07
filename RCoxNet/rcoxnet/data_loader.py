import os
import pandas as pd


def load_raw_data(root, cfg):
    data_dir = os.path.join(root, cfg['data_dir'])
    paths = {
        'mut':    os.path.join(data_dir, 'data_mutations.txt'),
        'clin_p': os.path.join(data_dir, 'data_clinical_patient.txt'),
        'clin_s': os.path.join(data_dir, 'data_clinical_sample.txt'),
        'ppi':    os.path.join(data_dir, cfg['ppi_file']),
    }

    print('Data files:')
    for label, path in paths.items():
        size = f'{os.path.getsize(path)/1e6:.1f} MB' if os.path.exists(path) else 'MISSING'
        print(f'  {label:8s}: {size}  ({path})')

    mut_raw    = pd.read_csv(paths['mut'],    sep='\t', low_memory=False)
    clin_p_raw = pd.read_csv(paths['clin_p'], sep='\t', skiprows=4)
    clin_s_raw = pd.read_csv(paths['clin_s'], sep='\t', skiprows=4)
    clin_p_raw.columns = clin_p_raw.columns.str.strip()
    clin_s_raw.columns = clin_s_raw.columns.str.strip()

    print(f'Mutations : {mut_raw.shape}  '
          f'patients: {mut_raw["Tumor_Sample_Barcode"].nunique()}')
    return mut_raw, clin_p_raw, clin_s_raw, paths['ppi'], None
