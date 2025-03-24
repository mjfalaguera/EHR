# EHR

## Objective

Connect UKBB and All of Us prescription data to ChEMBL and Open Targets Platform data (using IBD as a benchmark) in order to evaluate novel targets as therapeutic candidates.

## Data model

patient -> EHR -> prescription -> drug -> drug mechanism -> drug target VS novel target
![alt text](schema.png)

## Code

`queryData.py`: Query data from:
- ChEMBL: downloaded locally from https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_35/
- Open Targets Platform: downloaded locally from http://ftp.ebi.ac.uk/pub/databases/opentargets/platform/24.09
- Pathways: provided by Gowri Nayar
- UKBB

`nlpTools.py`: Apply NLP tools to drug data from UKBB and ChEMBL. Output stored in working directory `data/mix`.
