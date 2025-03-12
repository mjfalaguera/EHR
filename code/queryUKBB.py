#!/usr/bin/env python

__author__ = "Maria J. Falaguera"
__date__ = "12 Mar 2025"

"""
queryUKBB.py: Query EHR data for IBD patients in the UKBB (available at Stanford Sherlock cluster).
"""

import pandas as pd
import os

# EHR for IBD (Crohn's disease + Ulcerative colitis) patients
input1 = "/oak/stanford/groups/rbaltman/ukbiobank/metadata/coding19.tsv"
input2 = "/oak/stanford/groups/rbaltman/ukbiobank/phenotype_data/field_keys/extractedCols_everything.tsv"
output = "/home/groups/rbaltman/cote/data/ibd_extractedCols_everything.tsv"
if not os.path.exists(output):
    icd = pd.read_csv(input1, sep="\t")
    selectedCodes = (
        icd[
            icd.meaning.str.contains("K50")  # Crohn's disease"
            | icd.meaning.str.contains("K51")  # Ulcerative colitis
        ]
        .coding.unique()
        .tolist()
    )

    # get IBD participants
    with open(input2, "r") as fh:
        header = next(fh)
        eids = set()
        icd10fields = [
            idx
            for idx, field in enumerate(header.strip("\n").split("\t"))
            if "41202" in field  # icd fields prefix
        ]

        for line in fh:
            line = line.strip("\n").split("\t")
            for field in icd10fields:
                if line[field] in selectedCodes:
                    print(line[0])
                    eids.update([line[0]])
                    break

    # get EHR fields for IBD participants
    with open(input2, "r") as fhi, open(output, "w") as fho:
        header = next(fhi)
        fho.write(header)
        for line in fhi:
            if line.strip("\n").split("\t")[0] in eids:
                fho.write(line)

# EHR for IBD (Crohn's disease + Ulcerative colitis) patients
input1 = "/oak/stanford/groups/rbaltman/ukbiobank/phenotype_data/clinical_data/primary_care/gp_scripts.txt"
input2 = "/home/groups/rbaltman/cote/data/ibd_extractedCols_everything.tsv"
output = "/home/groups/rbaltman/cote/data/ibd_gp_scripts.tsv"
if not os.path.exists(output):
    with open(input1, "r", encoding="iso-8859-1") as fhi_gp_scripts, open(
        input2, "r"
    ) as fhi_ibd_fields, open(output, "w") as fho:

        header = next(fhi_ibd_fields)
        selectedParticipants = set(
            [line.split("\t")[0] for line in fhi_ibd_fields]
        )  # field 0: eid (participant)
        print(list(selectedParticipants)[:5])
        fho.write(next(fhi_gp_scripts))  # header

        for line in fhi_gp_scripts:
            participant = line.split("\t")[0]
            if participant in selectedParticipants:
                fho.write(line)

# IBD patient -> ICD
input = "/home/groups/rbaltman/cote/data/ibd_extractedCols_everything.tsv"
output = "/home/groups/rbaltman/cote/data/eidToIcd.tsv"  # move to working directory data/ukbb
if not os.path.exists(output):
    with open(input, "r") as fhi, open(output, "w") as fho:
        header = {
            col: idx
            for idx, col in enumerate(next(fhi).strip("\n").split("\t"))
            if "f.41202" in col
        }
        for line in fhi:
            line = line.strip("\n").split("\t")
            eid = line[0]
            icds = sorted(
                set([line[col] for col in header.values() if line[col] != "NA"])
            )
            for icd in icds:
                fho.write(
                    "{}\t{}\n".format(eid, icd)
                )  # one row per each patient -> icd

# IBD patient -> prescription
input = "/home/groups/rbaltman/cote/data/ibd_gp_scripts.tsv"
output = "/home/groups/rbaltman/cote/data/eidToPrescription.tsv"  # move to working directory data/ukbb
if not os.path.exists(output):
    data = pd.read_csv(input, sep="\t")[
        ["eid", "drug_name"]  # one row per each patient -> drug
    ].drop_duplicates()
    data.drug_name = data.apply(
        lambda row: row.drug_name.strip(), axis=1
    )  # some prescriptions have new line characters at the end -> remove them
    data.to_csv(output, sep="\t", index=False, header=False)
