#!/usr/bin/env python

__author__ = "Maria J. Falaguera"
__date__ = "18 Mar 2025"

"""
queryData.py: Query data from:
- ChEMBL: downloaded locally from https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_35/
- UKBB: downloaded locally from Stanford Sherlock cluster
- Open Targets Platform: downloaded locally from http://ftp.ebi.ac.uk/pub/databases/opentargets/platform/24.09
- Pathways: provided by Gowri Nayar

Output stored in working data/ directory.
"""

import os
import pandas as pd
import numpy as np
import oracledb
from pyspark.sql import functions as F
from pyspark.sql import SparkSession


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


def getInfoForMolecule(output=None, advanced=True):
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
    ]
    FROM = ["MOLECULE_DICTIONARY"]
    JOIN = [
        "MOLECULE_SYNONYMS ON MOLECULE_SYNONYMS.MOLREGNO = MOLECULE_DICTIONARY.MOLREGNO",
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


otPlatformDataDir = "/Users/mariaf/OT_platform/24.09/"
otDataDir = "../data/ot/"
chemblDataDir = "../data/chembl/"
mixDataDir = "../data/mix/"
ukbbDataDir = "../data/ukbb/"
pathwaysDataDir = "../data/pathways/"
reactomeDataDir = "../data/reactome/"

# IBD drug names in ChEMBL
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


# IBD drug mechanism of action in ChEMBL
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

# novel targets associated with IBD according to https://www.researchsquare.com/article/rs-5669559/v1
input1 = otPlatformDataDir + "evidenceIndirect/"
input2 = (
    "/Users/mariaf/TargetEngine/results/23.06/associationByDatasourceIndirectOverYears"
)
input3 = otPlatformDataDir + "targets"
input4 = otPlatformDataDir + "targetPrioritisation"
output = otDataDir + "novelTarget-EFO_0003767.tsv"
if not os.path.exists(output):
    with SparkSession.builder.getOrCreate() as spark:
        drugTarget = (
            spark.read.parquet(input1)
            .filter(F.col("diseaseId") == "EFO_0003767")
            .filter((F.col("datasourceId") == "chembl"))
            .withColumn("source", F.explode("urls.niceName"))
            # .filter(
            #     F.col("source").isin(
            #         ["DailyMed", "FDA", "EMA"]
            #     )  # only sources with approved drugs (clinical phase IV)
            # )
            # .select("targetId", "source", "drugId")
            .distinct()
        )

        (
            spark.read.parquet(input2)
            .filter(
                (F.col("diseaseId") == "EFO_0003767")
            )  # IBD (ontology propagation already applied)
            .filter(
                (F.col("novelty") > 0)
                & (F.col("novelty") >= 0.1)
                & (F.col("year") == 2023)
                & (F.col("datasourceId") != "europepmc")
                & (F.col("datasourceId") != "chembl")
            )
            .join(
                spark.read.parquet(input3)
                .withColumn("targetUniprot", F.explode("proteinIds"))
                .filter(F.col("targetUniprot.source") == "uniprot_swissprot")
                .select(
                    F.col("id").alias("targetId"),
                    F.col("approvedSymbol").alias("targetSymbol"),
                    F.col("approvedName").alias("targetName"),
                    F.col("biotype").alias("targetBiotype"),
                    F.col("targetUniprot.id").alias("targetUniprot"),
                ),
                on="targetId",
                how="left",
            )
            .join(drugTarget.select("targetId"), "targetId", "anti")
            # target prioritisation factors
            .join(spark.read.parquet(input4), how="left", on="targetId")
            .toPandas()
            .to_csv(output, sep="\t", index=False)
        )

# unique uniprots for novel targets
input = otDataDir + "novelTarget-EFO_0003767.tsv"
output = otDataDir + "novelTarget-EFO_0003767-targetUniprot.tsv"
if not os.path.exists(output):
    pd.read_csv(input, sep="\t")[["targetUniprot"]].drop_duplicates().to_csv(
        output, sep="\t", index=False, header=False
    )

# unique uniprots+ensembl for novel targets
input = otDataDir + "novelTarget-EFO_0003767.tsv"
output = otDataDir + "novelTarget-EFO_0003767-targetUniprot+targetId.tsv"
if not os.path.exists(output):
    pd.read_csv(input, sep="\t")[
        ["targetUniprot", "targetId"]
    ].drop_duplicates().to_csv(output, sep="\t", index=False, header=False)

# IBD and descedant diseases in OT
f = otDataDir + "disease.tsv"
if not os.path.exists(f):
    with SparkSession.builder.getOrCreate() as spark:
        spark.read.parquet("/Users/mariaf/OT_platform/24.09/diseases").filter(
            F.col("id") == "EFO_0003767"
        ).select(
            F.explode(F.array_union(F.col("descendants"), F.array(F.col("id")))).alias(
                "diseaseId"
            )
        ).join(
            spark.read.parquet("/Users/mariaf/OT_platform/24.09/diseases").select(
                F.col("id").alias("diseaseId"), F.col("name").alias("diseaseName")
            ),
            on="diseaseId",
            how="left",
        ).toPandas()[
            ["diseaseId", "diseaseName"]
        ].drop_duplicates().to_csv(
            f, sep="\t", index=False, header=False
        )

# IBD drugs in OT
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

# drug targets for IBD: IBD patient vs drug vs target
input1 = ukbbDataDir + "eidToIcd.tsv"
input2 = ukbbDataDir + "coding19.tsv"
input3 = ukbbDataDir + "eidToPrescription.tsv"
input4 = ukbbDataDir + "prescription.tsv"
input5 = ukbbDataDir + "drugName.tsv"
input6 = mixDataDir + "filteredClosestChemblDrugNameEmbForUkbbDrugNameEmb.tsv"
input7 = chemblDataDir + "drugId.tsv"
input10 = otPlatformDataDir + "targets"
input11 = otPlatformDataDir + "targetPrioritisation"
output = mixDataDir + "drugTarget-EFO_0003767.tsv"
if not os.path.exists(output):
    data = (
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
        # drug -> pref name
        .merge(
            getInfoForMolecule()[["CHEMBL_ID", "PREF_NAME"]].rename(
                columns={"CHEMBL_ID": "chemblId", "PREF_NAME": "chemblPrefName"}
            ),
            on="chemblId",
            how="left",
        )
        # drug -> target
        .merge(
            getTargetForDrug()[["CHEMBL_ID", "ACCESSION", "COMPONENT_SYNONYM"]]
            .dropna()
            .rename(
                columns={
                    "CHEMBL_ID": "chemblId",
                    "ACCESSION": "targetUniprot",
                    "COMPONENT_SYNONYM": "targetSymbol",
                }
            ),
            on="chemblId",
            how="left",
        )
        # drug -> mechanism
        .merge(
            getMechanismForDrug()[["CHEMBL_ID", "MECHANISM_OF_ACTION"]]
            .dropna()
            .rename(
                columns={"CHEMBL_ID": "chemblId", "MECHANISM_OF_ACTION": "mechanism"}
            ),
            on="chemblId",
            how="left",
        )
        .drop_duplicates()
    )

    with SparkSession.builder.getOrCreate() as spark:
        (
            # target info
            spark.read.parquet(input10)
            .select(
                F.col("id").alias("targetId"),
                F.col("approvedName").alias("targetName"),
                F.explode("proteinIds").alias("proteinId"),
            )
            .filter(F.col("proteinId.source") == "uniprot_swissprot")
            .select(
                "targetId",
                F.col("proteinId.id").alias("targetUniprot"),
                "targetName",
            )
            .filter(F.col("targetUniprot").isin(data.targetUniprot.unique().tolist()))
            # target prioritisation
            .join(spark.read.parquet(input11), on="targetId", how="left")
            .toPandas()
            .merge(
                data,
                on="targetUniprot",
                how="right",
            )
        ).to_csv(output, index=False, sep="\t")

# unique drug targets for IBD
input = mixDataDir + "drugTarget-EFO_0003767.tsv"
output = mixDataDir + "drugTarget-EFO_0003767-targetUniprot.tsv"
if not os.path.exists(output):
    pd.read_csv(input, sep="\t")[["targetUniprot"]].drop_duplicates().dropna().to_csv(
        output, sep="\t", header=False, index=False
    )

# pathways (modify  code as PPI)
input1 = mixDataDir + "drugTarget-EFO_0003767.tsv"
input2 = otDataDir + "novelTarget-EFO_0003767.tsv"
input3 = reactomeDataDir + "UniProt2Reactome_All_Levels.txt"
output = mixDataDir + "drugTarget-novelTarget-EFO_0003767-pathway.csv"
if not os.path.exists(output):

    pathway = (
        pd.read_csv(
            input3,
            sep="\t",
            names=["targetUniprot", "a", "b", "pathwayName", "c", "d"],
        )[["targetUniprot", "pathwayName"]]
        .dropna()
        .drop_duplicates()
    )

    drug = (
        pd.read_csv(input1, sep="\t")[["targetId", "targetUniprot", "targetSymbol"]]
        .dropna()
        .drop_duplicates()
        .assign(group="drug")
        .merge(pathway, on="targetUniprot", how="left")
    )

    novel = (
        pd.read_csv(input2, sep="\t")[["targetId", "targetUniprot", "targetSymbol"]]
        .dropna()
        .drop_duplicates()
        .assign(group="novel")
        .merge(pathway, on="targetUniprot", how="left")
    )

    pd.concat(
        [
            # drug-drug
            drug.rename(
                columns={c: c + "A" for c in drug.columns if c != "pathwayName"}
            ).merge(
                drug.rename(
                    columns={c: c + "B" for c in drug.columns if c != "pathwayName"}
                ),
                on="pathwayName",
                how="left",
            ),
            # novel-novel
            novel.rename(
                columns={c: c + "A" for c in novel.columns if c != "pathwayName"}
            ).merge(
                novel.rename(
                    columns={c: c + "B" for c in novel.columns if c != "pathwayName"}
                ),
                on="pathwayName",
                how="left",
            ),
            # drug-novel
            drug.rename(
                columns={c: c + "A" for c in drug.columns if c != "pathwayName"}
            ).merge(
                novel.rename(
                    columns={c: c + "B" for c in novel.columns if c != "pathwayName"}
                ),
                on="pathwayName",
                how="left",
            ),
            # novel-drug
            novel.rename(
                columns={c: c + "A" for c in novel.columns if c != "pathwayName"}
            ).merge(
                drug.rename(
                    columns={c: c + "B" for c in drug.columns if c != "pathwayName"}
                ),
                on="pathwayName",
                how="left",
            ),
        ]
    ).drop_duplicates().to_csv(
        output,
        sep="\t",
        index=False,
    )

# PPI (use evidence file instead to filter only direct ones)
input1 = mixDataDir + "drugTarget-EFO_0003767.tsv"
input2 = otDataDir + "novelTarget-EFO_0003767.tsv"
input3 = otPlatformDataDir + "interaction"
output = mixDataDir + "drugTarget-novelTarget-EFO_0003767-interaction.csv"
# if not os.path.exists(output):
if 1:
    with SparkSession.builder.getOrCreate() as spark:
        targets = spark.createDataFrame(
            pd.concat(
                [
                    # drug targets
                    pd.read_csv(input1, sep="\t")[
                        ["targetId", "targetUniprot", "targetSymbol"]
                    ]
                    .dropna()
                    .drop_duplicates()
                    .assign(group="drug"),
                    # novel targets
                    pd.read_csv(input2, sep="\t")[
                        ["targetId", "targetUniprot", "targetSymbol"]
                    ]
                    .dropna()
                    .drop_duplicates()
                    .assign(group="novel"),
                ]
            )
        )

        interaction = (
            spark.read.parquet(input3)
            .filter(F.col("sourceDatabase") != "reactome")
            .filter(F.col("scoring") == 1)
            .select(
                F.col("targetA").alias("targetIdA"),
                F.col("targetB").alias("targetIdB"),
                "sourceDatabase",
                "scoring",
            )
        )

        (
            targets.toDF(*[col + "A" for col in targets.columns])
            .join(interaction, on="targetIdA", how="left")
            .join(
                targets.toDF(*[col + "B" for col in targets.columns]),
                on="targetIdB",
                how="left",
            )
            .filter(F.col("targetIdA") != F.col("targetIdB"))
            # .withColumn("pair", F.array_sort(F.array(F.col("target"), col("last_name")))
            .distinct()
            .toPandas()
        ).to_csv(output, sep="\t", index=False)
