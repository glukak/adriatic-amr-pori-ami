# Please do not use this code before publication

# PORI and AMI: ecology-based indices of antimicrobial resistance hazard (Adriatic coast)

Code and data to reproduce the tables and figures of the manuscript (TO-FILL WHEN PUBLISHED: title, journal, DOI).

[!\[Binder](https://mybinder.org/badge\_logo.svg)](https://mybinder.org/v2/gh/USER/REPO/HEAD?labpath=run_pipeline.ipynb)
[!\[DOI](https://zenodo.org/badge/DOI/TO-FILL.svg)](https://doi.org/TO-FILL)

## Contents

* `adriatic\_amr\_pipeline.py`: complete pipeline (PORI, AMI, sentinel screen, robustness checks, figures).
* `data/adriaticlast.xlsx`: input data (36 stations, 17 watersheds, one campaign, 20-23 September 2024).
* `run\_pipeline.ipynb`: notebook that runs the pipeline, checks key numbers and shows the figures.
* `verifica.py`: compares key results with the manuscript.
* `requirements.txt`: package versions used.

## Run

Online: click the Binder badge, then *Run > Run All Cells* (about 1 minute; sessions are temporary, so download `results/`).

Locally:

```
python -m venv .venv \&\& source .venv/bin/activate
pip install -r requirements.txt
AMR\_XLSX=data/adriaticlast.xlsx AMR\_OUT=results python adriatic\_amr\_pipeline.py
python verifica.py
```

Outputs: `results/AMR\_results.xlsx` (all tables) and `results/figures/`. Random seed fixed (42).

## Licence

Code: MIT. Data: CC BY 4.0.

