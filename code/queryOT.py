#!/usr/bin/env python

__author__ = "Maria J. Falaguera"
__date__ = "12 Mar 2025"

"""
queryOT.py: Query Open Targets Platform data (downloaded locally from http://ftp.ebi.ac.uk/pub/databases/opentargets/platform/24.09). Output stored in working directory data/ot.
"""

from pyspark.sql import functions as F
from pyspark.sql import SparkSession
import os

dataDir = "../data/ot/"

spark = SparkSession.builder.getOrCreate()

# novel targets associated with IBD according to https://www.researchsquare.com/article/rs-5669559/v1
input1 = "/Users/mariaf/OT_platform/24.09/evidence/"
input2 = (
    "/Users/mariaf/TargetEngine/results/23.06/associationByDatasourceIndirectOverYears"
)
input3 = "/Users/mariaf/OT_platform/24.09/targets"
input4 = "/Users/mariaf/OT_platform/24.09/targetPrioritisation"
output = dataDir + "novelTarget.tsv"
if not os.path.exists(output):
    drugTarget = (
        spark.read.parquet(input1)
        .filter((F.col("datasourceId") == "chembl"))
        .withColumn("source", F.explode("urls.niceName"))
        .filter(
            F.col("source").isin(
                ["DailyMed", "FDA", "EMA"]
            )  # only sources with approved drugs (clinical phase IV)
        )
        .select("targetId", "source", "drugId")
        .distinct()
    )

    (
        spark.read.parquet(input2)
        .filter(
            (F.col("diseaseId") == "EFO_0003767")
        )  # IBD (ontology propagation already applied)
        .filter(
            (F.col("novelty") > 0)
            & (F.col("year") == 2023)
            & (F.col("datasourceId") != "europepmc")
        )
        .join(
            spark.read.parquet(input3)
            .withColumn("targetProteinId", F.explode("proteinIds"))
            .filter(F.col("targetProteinId.source") == "uniprot_swissprot")
            .select(
                F.col("id").alias("targetId"),
                F.col("approvedSymbol").alias("targetSymbol"),
                F.col("approvedName").alias("targetName"),
                F.col("biotype").alias("targetBiotype"),
                F.col("targetProteinId.id").alias("proteinId"),
            ),
            on="targetId",
            how="left",
        )
        .join(spark.read.parquet(input4), "targetId", "left")
        .join(drugTarget.select("targetId"), "targetId", "anti")
    ).toPandas().to_csv(output, sep="\t", index=False)

# IBD and descedant diseases
f = dataDir + "disease.tsv"
if not os.path.exists(f):
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

# targets metadata
f = dataDir + "target.tsv"
if not os.path.exists(f):
    (
        spark.read.parquet("/Users/mariaf/OT_platform/24.09/targets")
        .withColumn("targetId", F.explode("proteinIds"))
        .filter(F.col("targetId.source") == "uniprot_swissprot")
        .select(
            F.col("targetId.id").alias("targetId"),
            F.col("approvedSymbol").alias("targetSymbol"),
            F.col("approvedName").alias("targetName"),
        )
        .distinct()
        .toPandas()
        .to_csv(f, sep="\t", index=False, header=False)
    )
