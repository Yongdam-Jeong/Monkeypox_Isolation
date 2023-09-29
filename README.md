# Monkeypox_Isolation

Modelling the effectiveness of an isolation strategy for managing mpox outbreaks with variable infectiousness profiles

Yong Dam Jeong, Takara Nishiyama, Hyeongki Park, Masahiro Ishikane, Noriko Iwamoto, Kazuyuki Aihara, Koichi Watashi, Eline Op de Coul, William S Hart, Robin N Thompson, Norio Ohmagari, Jacco Wallinga, Shingo Iwami and Fuminari Miura



1. Viral clearance model

'mpox.csv' is an example data of longitudinal viral load data for mpox patients.

'model.txt' is the description for mathematical model of viral clearance.

These files are used to estimate parameters of the viral clearance model in MONOLIX2023R1.



2. Simulation for isolation guideline

'populationParameters_mpox.txt' is an example for the estimated parameters of viral clearance model for mpox patients.

'Fixed-duration.R' calculates numerical values of "Risk of prematurely ending isolation", "Infectious period after ending isolation", and "Unnecessarily prolonged isolation period" for mpox patients under fixed-duration isolation rules.

'Testing-based.R' calculates numerical values of "Risk of prematurely ending isolation", "Infectious period after ending isolation", and "Unnecessarily prolonged isolation period" for mpox patients under testing-based isolation rules.

* Text files and R codes above should be in the same location.
