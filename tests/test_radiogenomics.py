import numpy as np
import pandas as pd

from mumoas.radiogenomics.association import fdr_bh, spearman_program_associations


def test_fdr_bh_returns_monotonic_adjusted_p_values():
    adjusted = fdr_bh([0.001, 0.02, 0.03, 0.8])

    assert np.all(np.diff(adjusted) >= 0)


def test_spearman_program_associations_returns_long_table_for_features_and_programs():
    imaging = pd.DataFrame(
        {
            "patient_id": ["P001", "P002", "P003", "P004"],
            "feature_a": [1.0, 2.0, 3.0, 4.0],
            "feature_b": [4.0, 3.0, 2.0, 1.0],
            "batch": ["A", "A", "B", "B"],
        }
    )
    programs = pd.DataFrame(
        {
            "patient_id": ["P001", "P002", "P003", "P004"],
            "MP1": [1.0, 2.0, 3.0, 4.0],
            "MP2": [4.0, 3.0, 2.0, 1.0],
        }
    )

    associations = spearman_program_associations(
        imaging=imaging,
        programs=programs,
        patient_col="patient_id",
        program_cols=["MP1", "MP2"],
    )

    assert list(associations.columns) == ["feature", "program", "rho", "p_value", "fdr"]
    assert len(associations) == 4
    assert set(associations["feature"]) == {"feature_a", "feature_b"}
    assert set(associations["program"]) == {"MP1", "MP2"}
