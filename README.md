# Mpox_Isolation

# Viral clearance model

'mpox.csv' is an example data of longitudinal viral load data for mpox patients.

'model.txt' is the description for mathematical model of viral clearance.

These files are used to estimate parameters of the viral clearance model in MONOLIX2023R1.


# Simulation for isolation guideline

'populationParameters_mpox.txt' is an example for the estimated parameters of viral clearance model for mpox patients.

'Symptom-based.R' calculates numerical values of "Risk of prematurely ending isolation", "Infectious period after ending isolation", and "Unnecessarily prolonged isolation period" for mpox patients under symptom-based isolation rules.

'Fixed-duration.R' calculates numerical values of "Risk of prematurely ending isolation", "Infectious period after ending isolation", and "Unnecessarily prolonged isolation period" for mpox patients under fixed-duration isolation rules.

'Testing-based.R' calculates numerical values of "Risk of prematurely ending isolation", "Infectious period after ending isolation", and "Unnecessarily prolonged isolation period" for mpox patients under testing-based isolation rules.

* Text files and R (version 4.2.3) codes above should be in the same location.
