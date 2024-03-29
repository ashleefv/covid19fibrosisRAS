# covid19fibrosisRAS
[![DOI](https://zenodo.org/badge/557388588.svg)](https://zenodo.org/doi/10.5281/zenodo.7295153)

## Overview 
The mathematical model describes the role of patient-specific premorbidity, age, and sex differences on Renin-Angiotensin System.

## Mathematical Modeling of Impacts of Patient Differences on Renin-Angiotensin System and Applications to COVID-19 Lung Fibrosis Outcomes
### Manuscript: [DOI 10.1101/2022.11.06.515367v2](https://doi.org/10.1101/2022.11.06.515367v2)
### Code Authors
Mohammad Aminul Islam and Ashlee N. Ford Versypt, 
Dept. of Chemical & Biological Engineering,
University at Buffalo, The State University of New York.

Corresponding author: A. N. Ford Versypt, ashleefv@buffalo.edu

### Scripts

* [RAS.ipynb] This file will generate RAS and immune dynamics and dose response to ueACE2

* [RAS_initial_value_variations.ipynb] This file generate prediction with variation of initial values of ANGI, ANGII, and ANG1-7 within experimentally obserevd SARS-CoV-2 negative and SARS-CoV-2 positive patients

* [RAS_local_sensitivity_parameters.ipynb] This file generate time dynamic sensitivity for model parameters

* [RAS_ANGII_to_TGF.ipynb] This file generate fibrosis model dynamics with variation of ANGII induced TGF-beta production rate (kT) for ueACE_1000

* ueACE2_200.p - ueACE2_2000.p and ueACE2t are data file for ACE2 dynamics from the agent-based model

* The folder [Template of in silico experiments] contains template code for all the cases of ACE variations from 200-2000

* [Multi_run_analysis.ipynb] This file extracts data from multiple replications and plots mean responses. This file also generates the ueACE2 dynamics


[RAS.ipynb]: https://github.com/ashleefv/covid19fibrosisRAS/blob/master/RAS.ipynb
[RAS_local_sensitivity_parameters.ipynb]: https://github.com/ashleefv/covid19fibrosisRAS/blob/master/RAS_local_sensitivity_parameters.ipynb
[RAS_initial_value_variations.ipynb]: https://github.com/ashleefv/covid19fibrosisRAS/blob/master/RAS_initial_value_variations.ipynb
[RAS_ANGII_to_TGF.ipynb]: https://github.com/ashleefv/covid19fibrosisRAS/blob/master/RAS_ANGII_to_TGF.ipynb
[Template of in silico experiments]: https://github.com/ashleefv/covid19fibrosisRAS/tree/master/Template%20of%20in%20silico%20experiments
[Multi_run_analysis.ipynb]: https://github.com/ashleefv/covid19fibrosisRAS/blob/master/Multi_run_analysis.ipynb


## Acknowledgements
Research reported in this publication was supported by the National Institute of General Medical Sciences of the National Institutes of Health under award number R35GM133763. The content is solely the responsibility of the authors and does not necessarily represent the offcial views of the National Institutes of Health.
