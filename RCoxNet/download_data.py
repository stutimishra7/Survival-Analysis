"""
Download TCGA clinical and mutation data from cBioPortal.

Usage:
    python download_data.py --cancer LUNG
    python download_data.py --cancer ALL     # LUNG + GBM + OV (BRCA already included)
"""

import argparse
import os
import sys
import urllib.request

BASE_URL = 'https://cbioportal-datahub.s3.amazonaws.com'

STUDIES = {
    'LUNG': 'luad_tcga_pan_can_atlas_2018',
    'GBM':  'gbm_tcga_pan_can_atlas_2018',
    'OV':   'ov_tcga_pan_can_atlas_2018',
}

FILES = [
    'data_mutations.txt',
    'data_clinical_patient.txt',
    'data_clinical_sample.txt',
]

ROOT = os.path.dirname(os.path.abspath(__file__))

parser = argparse.ArgumentParser()
parser.add_argument('--cancer', default='ALL',
                    choices=list(STUDIES) + ['ALL'])
args = parser.parse_args()

targets = list(STUDIES) if args.cancer == 'ALL' else [args.cancer]


def _progress(count, block_size, total):
    done = count * block_size
    if total > 0:
        pct = min(done / total * 100, 100)
        bar = '#' * int(pct / 2)
        sys.stdout.write(f'\r  [{bar:<50}] {pct:.0f}%')
        sys.stdout.flush()


def download_file(url, dest):
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    if os.path.exists(dest):
        print(f'  already exists: {os.path.basename(dest)}')
        return
    print(f'  {os.path.basename(dest)}')
    try:
        urllib.request.urlretrieve(url, dest, reporthook=_progress)
        print()
    except Exception as e:
        print(f'\n  failed: {e}')
        if os.path.exists(dest):
            os.remove(dest)


for cancer in targets:
    study_id = STUDIES[cancer]
    print(f'\n{cancer}  ({study_id})')
    for fname in FILES:
        url  = f'{BASE_URL}/{study_id}/{fname}'
        dest = os.path.join(ROOT, 'data', cancer, fname)
        download_file(url, dest)

print('\nDone. Place tissue-specific PPI files under data/<CANCER>/ — see data/README.md.')
