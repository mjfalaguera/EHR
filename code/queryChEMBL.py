#!/usr/bin/env python

__author__ = "Maria J. Falaguera"
__date__ = "12 Mar 2025"

"""
queryChEMBL.py: Query ChEMBL data (downloaded locally from https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_35/). Output stored in working directory data/chembl.
"""

import os
import pandas as pd
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


def getInfoForMolecule(output=None, advanced=False):
    """
    Get metadata for molecule.
    """
    SELECT = [
        "MOLECULE_DICTIONARY.CHEMBL_ID",
        "MOLECULE_DICTIONARY.FIRST_APPROVAL",
        "MOLECULE_DICTIONARY.PREF_NAME",
        "MOLECULE_DICTIONARY.MOLREGNO",
        "MOLECULE_DICTIONARY.MOLECULE_TYPE",
        "MOLECULE_SYNONYMS.SYNONYMS",
        "MOLECULE_ATC_CLASSIFICATION.LEVEL5 AS ATC",
    ]
    FROM = ["MOLECULE_DICTIONARY"]
    JOIN = [
        "MOLECULE_SYNONYMS ON MOLECULE_SYNONYMS.MOLREGNO = MOLECULE_DICTIONARY.MOLREGNO",
        "MOLECULE_ATC_CLASSIFICATION ON MOLECULE_ATC_CLASSIFICATION.MOLREGNO = MOLECULE_DICTIONARY.MOLREGNO",
    ]
    WHERE = ["1 = 1"]
    if advanced == True:
        WHERE.append("(MOLECULE_DICTIONARY.MAX_PHASE >= 3)")

    # build query
    query = "SELECT {} FROM {} LEFT JOIN {} WHERE {}".format(
        ", ".join(SELECT),
        ", ".join(FROM),
        " LEFT JOIN ".join(JOIN),
        " AND ".join(WHERE),
    )

    # parse results
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


def getTargetForDrug():
    """
    Get drug targets.
    """

    SELECT = [
        "MOLECULE_DICTIONARY.CHEMBL_ID",
        "COMPONENT_SEQUENCES.ACCESSION",
        "COMPONENT_SYNONYMS.COMPONENT_SYNONYM",
    ]
    FROM = ["DRUG_MECHANISM"]
    JOIN = [
        # mechanism propagated from children to parent to avoid parents with missing mechanism annotation
        "MOLECULE_HIERARCHY ON DRUG_MECHANISM.MOLREGNO=MOLECULE_HIERARCHY.MOLREGNO",
        "MOLECULE_DICTIONARY ON MOLECULE_DICTIONARY.MOLREGNO=MOLECULE_HIERARCHY.PARENT_MOLREGNO",
        "TARGET_COMPONENTS ON TARGET_COMPONENTS.TID=DRUG_MECHANISM.TID",
        "COMPONENT_SEQUENCES ON TARGET_COMPONENTS.COMPONENT_ID=COMPONENT_SEQUENCES.COMPONENT_ID",
        "COMPONENT_SYNONYMS ON COMPONENT_SYNONYMS.COMPONENT_ID=COMPONENT_SEQUENCES.COMPONENT_ID",
    ]
    WHERE = ["COMPONENT_SYNONYMS.SYN_TYPE = 'GENE_SYMBOL'"]

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


def getMechanismForDrug():
    """
    Get drug mechanism of action.
    """

    SELECT = [
        "MOLECULE_DICTIONARY.CHEMBL_ID",
        "DRUG_MECHANISM.MECHANISM_OF_ACTION",
    ]
    FROM = ["DRUG_MECHANISM"]
    JOIN = [
        # mechanism propagated from children to parent to avoid parents with missing mechanism annotation
        "MOLECULE_HIERARCHY ON DRUG_MECHANISM.MOLREGNO=MOLECULE_HIERARCHY.MOLREGNO",
        "MOLECULE_DICTIONARY ON MOLECULE_DICTIONARY.MOLREGNO=MOLECULE_HIERARCHY.PARENT_MOLREGNO",
    ]

    query = "SELECT {} FROM {} LEFT JOIN {}".format(
        ", ".join(SELECT),
        ", ".join(FROM),
        " LEFT JOIN ".join(JOIN),
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


chemblDataDir = "../data/chembl/"
mixDataDir = "../data/mix/"

# IBD drug names
input = mixDataDir + "diseaseIdToDrugId.tsv"
output1 = chemblDataDir + "drugId.tsv"
output2 = chemblDataDir + "drugName.tsv"
if (not os.path.exists(output1)) or (not os.path.exists(output2)):
    data = (
        getInfoForMolecule(advanced=True)[["CHEMBL_ID", "PREF_NAME", "SYNONYMS"]]
        .dropna()
        .rename(columns={"CHEMBL_ID": "drugId"})
        .merge(
            pd.read_csv(input, names=["diseaseId", "drugId"], sep="\t"),
            on="drugId",
            how="inner",
        )
    )

    data = pd.concat(
        [
            data.rename(columns={"PREF_NAME": "drugName"})[["drugId", "drugName"]],
            data.rename(columns={"SYNONYMS": "drugName"})[["drugId", "drugName"]],
        ]
    ).drop_duplicates()
    data[["drugId"]].to_csv(output1, sep="\t", index=False, header=False)
    data[["drugName"]].to_csv(output2, sep="\t", index=False, header=False)

# IBD drug targets
input = chemblDataDir + "drugId.tsv"
output = chemblDataDir + "drugIdToTargetId.tsv"
if not os.path.exists(output):
    (
        getTargetForDrug()[["CHEMBL_ID", "ACCESSION", "COMPONENT_SYNONYM"]]
        .dropna()
        .rename(
            columns={
                "CHEMBL_ID": "drugId",
                "ACCESSION": "targetId",
                "COMPONENT_SYNONYM": "targetSymbol",
            }
        )
        .merge(
            pd.read_csv(input, names=["drugId"], sep="\t"),
            on="drugId",
            how="inner",
        )
    ).drop_duplicates().to_csv(output, sep="\t", index=False, header=False)

# IBD drug mechanism of action
input = chemblDataDir + "drugId.tsv"
output = chemblDataDir + "drugIdToMechanism.tsv"
if not os.path.exists(output):
    (
        getMechanismForDrug()[["CHEMBL_ID", "MECHANISM_OF_ACTION"]]
        .dropna()
        .rename(columns={"CHEMBL_ID": "drugId", "MECHANISM_OF_ACTION": "mechanism"})
        .merge(
            pd.read_csv(input, names=["drugId"], sep="\t"),
            on="drugId",
            how="inner",
        )
    ).drop_duplicates().to_csv(output, sep="\t", index=False, header=False)
