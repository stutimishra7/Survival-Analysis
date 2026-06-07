# Data

BRCA data is included in this repository (mutations via Git LFS).
For LUNG, GBM, and OV, download the files from cBioPortal as described below.

Expected layout:

```
data/
  BRCA/
    data_mutations.txt          (included, via Git LFS)
    data_clinical_patient.txt   (included)
    data_clinical_sample.txt    (included)
    breast-PPI.txt              (included)
  LUNG/
    data_mutations.txt
    data_clinical_patient.txt
    data_clinical_sample.txt
    lung-PPI.txt
  GBM/
    data_mutations.txt
    data_clinical_patient.txt
    data_clinical_sample.txt
    glioblastoma-PPI.txt
  OV/
    data_mutations.txt
    data_clinical_patient.txt
    data_clinical_sample.txt
    ovarian-PPI-RALF.txt
```

---

## Downloading LUNG / GBM / OV data

Run the provided script:

```bash
python download_data.py --cancer LUNG
python download_data.py --cancer GBM
python download_data.py --cancer OV
```

This fetches the three clinical/mutation files per cancer from cBioPortal
(TCGA PanCancer Atlas 2018 studies).

---

## PPI files

Tissue-specific PPI files were derived from ConsensusPathDB v35
(https://cpdb.molgen.mpg.de/). The `breast-PPI.txt` for BRCA is included.

For LUNG, GBM, and OV PPI files, contact the authors:
stutik@iiitd.ac.in
