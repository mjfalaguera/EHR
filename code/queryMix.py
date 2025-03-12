#!/usr/bin/env python

__author__ = "Maria J. Falaguera"
__date__ = "12 Mar 2025"

"""
queryMix.py: Query data from ChEMBL, UKBB, Open Targets Platform, pathways. Output stored in working directory data/mix.
"""

import os
import pandas as pd
import sys
import numpy as np
import oracledb


def openConnection():
    """
    Connect to ChEMBL database.
    """
    service = "chempro"
    user = "chembl_34"
    host = "ora-vm-089.ebi.ac.uk"
    password = "chembl_34"
    port = 1531

    dsn = f"{host}:{port}/{service}"
    conn = oracledb.connect(user=user, password=password, dsn=dsn)

    return conn


def getIndicationForDrug(propagation=True):
    """
    Get drug indication.
    """
    SELECT = [
        "MOLECULE_DICTIONARY.CHEMBL_ID",
        "DRUG_INDICATION.EFO_ID",
    ]
    FROM = [
        "MOLECULE_HIERARCHY"
    ]  # mechanism propagated from children to parent to avoid parents with missing mechanism annotation
    JOIN = [
        "DRUG_INDICATION ON DRUG_INDICATION.MOLREGNO=MOLECULE_HIERARCHY.MOLREGNO",
        "MOLECULE_DICTIONARY ON MOLECULE_DICTIONARY.MOLREGNO=MOLECULE_HIERARCHY.PARENT_MOLREGNO",
    ]
    WHERE = ["DRUG_INDICATION.MAX_PHASE_FOR_IND >= 3"]

    query = "SELECT {} FROM {} LEFT JOIN {} WHERE {}".format(
        ", ".join(SELECT),
        ", ".join(FROM),
        " LEFT JOIN ".join(JOIN),
        " AND ".join(WHERE),
    )

    results = []
    with openConnection() as conn:
        with conn.cursor() as cur:
            cur.execute(query)
            cur.rowfactory = lambda *args: dict(
                zip([d[0] for d in cur.description], args)
            )
            res = cur.fetchall()
            for r in res:
                results.append(r)
    results = pd.DataFrame(results)
    return results


otDataDir = "../data/ot/"
chemblDataDir = "../data/chembl/"
mixDataDir = "../data/mix/"
ukbbDataDir = "../data/ukbb/"
pathwaysDataDir = "../data/pathways/"

# IBD drugs
input = otDataDir + "disease.tsv"
output = mixDataDir + "diseaseIdToDrugId.tsv"
if not os.path.exists(output):
    data = (
        getIndicationForDrug()
        .rename(columns={"EFO_ID": "diseaseId", "CHEMBL_ID": "chemblId"})[
            ["diseaseId", "chemblId"]
        ]
        .drop_duplicates()
    )
    data.diseaseId = data.diseaseId.str.replace(":", "_")
    pd.read_csv(input, sep="\t", names=["diseaseId"]).merge(
        data,
        on="diseaseId",
        how="inner",
    ).to_csv(output, index=False, header=False, sep="\t")

# closest match in ChEMBL found for UKBB drug name
input1 = mixDataDir + "closestChemblDrugNameEmbForUkbbDrugNameEmb.npy"
input2 = chemblDataDir + "drugName.tsv"
input3 = ukbbDataDir + "drugName.tsv"
output = mixDataDir + "niceClosestChemblDrugNameEmbForUkbbDrugNameEmb.tsv"
if not os.path.exists(output):
    (
        pd.DataFrame(np.load(input1))
        .reset_index()
        .rename(columns={"index": "ukbbDrugId", 0: "chemblDrugId", 1: "distance"})
        .merge(
            pd.read_csv(input2, sep="\t", names=["drugName"])
            .reset_index()
            .rename(columns={"index": "chemblDrugId"}),
            on="chemblDrugId",
            how="left",
        )
        .merge(
            pd.read_csv(input3, sep="\t", names=["drugName"])
            .reset_index()
            .rename(columns={"index": "ukbbDrugId"}),
            on="ukbbDrugId",
            how="left",
        )
    ).sort_values(by="distance").to_csv(output, sep="\t", header=False, index=False)

# dataframe of closest match in ChEMBL found for UKBB drug name
input = mixDataDir + "niceClosestChemblDrugNameEmbForUkbbDrugNameEmb.tsv"
output = mixDataDir + "filteredClosestChemblDrugNameEmbForUkbbDrugNameEmb.tsv"
if not os.path.exists(output):
    data = pd.read_csv(
        input,
        sep="\t",
        names=[
            "ukbbDrugId",
            "chemblDrugId",
            "distance",
            "chemblDrugName",
            "ukbbDrugName",
        ],
    )
    data[data.distance < 0.001].to_csv(output, sep="\t", header=False, index=False)

# triads: IBD patient vs drug vs target
input1 = ukbbDataDir + "eidToIcd.tsv"
input2 = ukbbDataDir + "coding19.tsv"
input3 = ukbbDataDir + "eidToPrescription.tsv"
input4 = ukbbDataDir + "prescription.tsv"
input5 = ukbbDataDir + "drugName.tsv"
input6 = mixDataDir + "filteredClosestChemblDrugNameEmbForUkbbDrugNameEmb.tsv"
input7 = chemblDataDir + "drugId.tsv"
input8 = chemblDataDir + "drugIdToTargetId.tsv"
input9 = chemblDataDir + "drugIdToMechanism.tsv"
output = mixDataDir + "triads.tsv"
if not os.path.exists(output):
    (
        # eid -> icd
        pd.read_csv(
            input1,
            sep="\t",
            names=[
                "eid",
                "icd",
            ],
        )
        # icd -> meaning
        .merge(
            pd.read_csv(input2, sep="\t").rename(columns={"coding": "icd"})[
                ["icd", "meaning"]
            ],
            on="icd",
            how="left",
        )[["eid", "meaning"]]
        .rename(columns={"meaning": "icd"})
        # eid -> prescription
        .merge(
            pd.read_csv(input3, sep="\t", names=["eid", "prescription"]),
            on="eid",
            how="left",
        )
        # prescription -> drugName
        .merge(
            pd.read_csv(input4, sep="\t", names=["prescription"])
            .reset_index()
            .rename(columns={"index": "ukbbDrugId"}),
            on="prescription",
            how="left",
        )
        .merge(
            pd.read_csv(input5, sep="\t", names=["ukbbDrugName"])
            .reset_index()
            .rename(columns={"index": "ukbbDrugId"}),
            on="ukbbDrugId",
            how="left",
        )
        # prescription -> chemblId
        .merge(
            pd.read_csv(
                input6,
                sep="\t",
                names=[
                    "ukbbDrugId",
                    "chemblDrugId",
                    "distance",
                    "chemblDrugName",
                    "ukbbDrugName",
                ],
            ),
            on=["ukbbDrugId", "ukbbDrugName"],
            how="left",
        )
        .merge(
            pd.read_csv(
                input7,
                sep="\t",
                names=[
                    "chemblId",
                ],
            )
            .reset_index()
            .rename(columns={"index": "chemblDrugId"}),
            on="chemblDrugId",
            how="left",
        )
        # drug -> target
        .merge(
            pd.read_csv(
                input8, sep="\t", names=["chemblId", "targetId", "targetSymbol"]
            ),
            on="chemblId",
            how="left",
        )
        # drug -> mechanism
        .merge(
            pd.read_csv(input9, sep="\t", names=["chemblId", "mechanism"]),
            on="chemblId",
            how="left",
        )
    ).to_csv(output, index=False, sep="\t")

# work in progress:
# input = pathwaysDataDir + "novel_target_dist.csv"
# output = mixDataDir + "triads.tsv"
# if not os.path.exists(output):
#     pd.read_csv(input, index_col=0).reset_index().rename(columns={"index": "targetId"})
