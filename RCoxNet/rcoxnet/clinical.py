"""
Clinical feature extraction for TCGA cohorts.

Merges patient-level survival data with sample-level TMB and MSI scores.
GBM patients with missing age are excluded per the paper's preprocessing.
"""

import numpy as np
import pandas as pd

from rcoxnet.config import CLIN_FEATURES


def extract_clinical(clin_p_raw, clin_s_raw, seed_dict, mut_count_dict, cfg):
    # patient-level file: survival endpoints and age
    clin_p = clin_p_raw.copy()
    clin_p['OS_STATUS'] = (clin_p['OS_STATUS'].astype(str)
                           .str.split(':').str[0].str.strip()
                           .replace('', np.nan).astype(float))
    clin_p['OS_MONTHS'] = pd.to_numeric(clin_p['OS_MONTHS'], errors='coerce')
    clin_p['AGE']       = pd.to_numeric(clin_p['AGE'],       errors='coerce')
    patient_clin = clin_p[['PATIENT_ID', 'OS_STATUS', 'OS_MONTHS', 'AGE']].copy()

    # sample-level file: TMB and MSI scores
    clin_s = clin_s_raw.copy()
    clin_s['tmb_score'] = pd.to_numeric(
        clin_s.get('TMB_NONSYNONYMOUS', np.nan), errors='coerce')
    clin_s['msi_score'] = pd.to_numeric(
        clin_s.get(cfg['msi_col'], np.nan), errors='coerce')
    # keep primary tumour samples only (barcode ends with -01)
    clin_s = clin_s[clin_s['SAMPLE_ID'].astype(str).str.endswith('-01')].copy()
    sample_clin = clin_s[['PATIENT_ID', 'SAMPLE_ID', 'tmb_score', 'msi_score']].copy()

    # merge and drop patients missing survival information
    clinical = patient_clin.merge(sample_clin, on='PATIENT_ID', how='inner')
    clinical = (clinical
                .dropna(subset=['OS_MONTHS', 'OS_STATUS'])
                .reset_index(drop=True))

    if cfg['drop_missing_age']:
        before   = len(clinical)
        clinical = clinical.dropna(subset=['AGE']).reset_index(drop=True)
        print(f'Dropped {before - len(clinical)} patients with missing AGE')

    # TCGA barcodes are sample-level (TCGA-XX-XXXX-01); trim to patient level
    seed_by_patient = {}
    for barcode, genes in seed_dict.items():
        pid = '-'.join(str(barcode).split('-')[:3])
        seed_by_patient.setdefault(pid, set()).update(genes)
    seed_by_patient = {k: sorted(v) for k, v in seed_by_patient.items()}

    mut_count_by_patient = {}
    for barcode, gene_counts in mut_count_dict.items():
        if not isinstance(gene_counts, dict):
            continue
        pid = '-'.join(str(barcode).split('-')[:3])
        for g, c in gene_counts.items():
            mut_count_by_patient.setdefault(pid, {})
            mut_count_by_patient[pid][g] = mut_count_by_patient[pid].get(g, 0) + c

    both     = set(clinical['PATIENT_ID']) & set(seed_by_patient)
    clinical = clinical[clinical['PATIENT_ID'].isin(both)].reset_index(drop=True)

    print(f'Final cohort     : {len(clinical)} patients')
    print(f'Event rate       : {clinical["OS_STATUS"].mean():.2%}')
    print(f'Clinical features: {CLIN_FEATURES}')
    return clinical, seed_by_patient, mut_count_by_patient
