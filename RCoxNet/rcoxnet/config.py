# clinical features passed to the model — same for all cancer types
CLIN_FEATURES = ['AGE', 'msi_score', 'tmb_score']

# columns kept in the merged feature dataframe
META_COLS = ['SAMPLE_ID', 'OS_MONTHS', 'OS_STATUS', 'AGE', 'tmb_score', 'msi_score']

# per-cancer data paths and preprocessing flags
CANCER_CONFIGS = {
    'LUNG': {
        'data_dir':         'data/LUNG',
        'ppi_file':         'lung-PPI.txt',
        'results_dir':      'Results/Pipeline_LUNG',
        'msi_col':          'MSI_SENSOR_SCORE',
        'drop_missing_age': False,
    },
    'OV': {
        'data_dir':         'data/OV',
        'ppi_file':         'ovarian-PPI-RALF.txt',
        'results_dir':      'Results/Pipeline_OV',
        'msi_col':          'MSI_SENSOR_SCORE',
        'drop_missing_age': False,
    },
    'GBM': {
        'data_dir':         'data/GBM',
        'ppi_file':         'glioblastoma-PPI.txt',
        'results_dir':      'Results/Pipeline_GBM',
        'msi_col':          'MSI_SCORE_MANTIS',
        'drop_missing_age': True,
    },
    'BRCA': {
        'data_dir':         'data/BRCA',
        'ppi_file':         'breast-PPI.txt',
        'results_dir':      'Results/Pipeline_BRCA',
        'msi_col':          'MSI_SENSOR_SCORE',
        'drop_missing_age': False,
    },
}

# default training settings
N_REPEATS    = 20
TEST_FRAC    = 0.20
VAL_FRAC     = 0.20
EPOCHS       = 1000
CHECKPOINT   = 50
GRAD_CLIP    = 100

# hyperparameter search grid
GRID_ARCHITECTURES = [(8, 4, 1), (32, 16, 8), (128, 64, 16), (64, 16, 4)]
GRID_L2            = [0.001, 0.01, 0.03, 0.1]
GRID_LR            = [1e-3, 3e-4, 1e-4]
GS_SEED            = 42
GS_EPOCHS          = 600
GS_CHECKPOINT      = 50
GS_GRAD_CLIP       = 100

# RWR and feature selection settings
ALPHA_VALUES = [0.1, 0.3, 0.5, 0.7, 0.8]
RWR_MAX      = 200
RWR_TOL      = 1e-10
