HPO analyses
---
1. Use `snapshot_patient_hpo.py` to take a snapshot of the patient database with hpos. This will avoid errors during analyses when patients' data are modified at client end. Note that each patient's hpo terms are expanded to include all ancestor terms [with no repeat].

2. Use `get_hpo_freq.py` to get hpo's frequencies. 

1 and 2 are also both needed for phenogenon.

3. Use `hpo_matrix.py` to produce a hpo co-occurrence matrix. 

4. Use `minimise_snapshot.py` to minise each patient's hpo terms. this is for downstream simulation

5. Use `simulate.py` 

