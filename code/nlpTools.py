#!/usr/bin/env python

__author__ = "Maria J. Falaguera"
__date__ = "12 Mar 2025"

"""
nlpTools.py: Apply NLP tools to drug data from UKBB and ChEMBL. Output stored in working directory data/mix.
"""

import os
import numpy as np
import pandas as pd
from transformers import (
    AutoModelForTokenClassification,
    AutoModel,
    AutoTokenizer,
    pipeline,
    GPT2LMHeadModel,
    GPT2Tokenizer,
)
import torch


def nerDrug(prescription, model="jsylee/scibert_scivocab_uncased-finetuned-ner"):
    """
    Extract drug names from prescription text.
    """
    if not isinstance(prescription, list):
        prescription = [prescription]

    prescription = [
        s.strip() + ";" for s in prescription
    ]  # I found that adding a delimiter after one-word drug names improves the model accuracy to recognise drug names // Strip to avoid inconsistencies later

    device = (
        torch.device(
            "mps"
        )  # MPS (Metal Performance Shaders) device, which is Apple's GPU framework for M1 Macs
        if torch.backends.mps.is_available()  # if MPS device is available in your mac and you don't move the model, inputs and outputs to it you get an error
        else torch.device("cpu")  # by default model is stored in CPU
    )
    tokenizer = AutoTokenizer.from_pretrained(
        model, truncation=True, model_max_length=512
    )
    model = AutoModelForTokenClassification.from_pretrained(
        model,
        num_labels=5,
        id2label={0: "O", 1: "B-DRUG", 2: "I-DRUG", 3: "B-EFFECT", 4: "I-EFFECT"},
    ).to(
        device
    )  # by default model is stored in CPU so we need to move it to GPU

    # NER
    modelPipeline = pipeline(
        task="ner",
        model=model,
        tokenizer=tokenizer,
        device=device,
    )
    ner = modelPipeline(prescription)

    # concat NER results
    drugResult = []
    for s in ner:
        drugs = set()
        drug = ""
        for entity in s:
            token = entity["word"]
            tokenType = entity["entity"]
            tokenScore = entity["score"]
            if (tokenScore >= 0.5) and ("DRUG" in tokenType):
                if (
                    "B-DRUG" in tokenType
                ):  # begining of drug (B-DRUG) or effect (B-EFFECT) -> I found some drug name classified as EFFECT
                    if (
                        drug != ""
                    ):  # at the beginning of a new drug name, if there was a previous drug already, save the previous one before processing the new one
                        drugs.update([drug])
                        drug = ""  # reinitiate drug name
                    drug += token
                elif ("I-DRUG" in tokenType) and (drug != ""):
                    if "##" in token:
                        token = token.strip("##")
                    else:
                        token = " " + token
                    drug += token
        drugs.update([drug])  # don't forget adding the last drug
        drugResult.append(drugs)

    # different drugs can be found for the same prescription
    return {
        "prescription": [
            p.strip(";") for p, drugs in zip(prescription, drugResult) for d in drugs
        ],
        "drug": [d for p, drugs in zip(prescription, drugResult) for d in drugs],
    }


def getEmbeddingForText(
    text,
    # model="emilyalsentzer/Bio_ClinicalBERT",
    model="dmis-lab/biobert-v1.1",
    batch_size=1000,
):
    """
    Calculate embedding for text.
    """
    if not isinstance(text, list):
        text = [text]

    text = [t.lower() for t in text]  # BERT is case-sensitive by default

    torch.mps.empty_cache()
    device = torch.device("mps" if torch.backends.mps.is_available() else "cpu")

    tokenizer = AutoTokenizer.from_pretrained(model)
    model = AutoModel.from_pretrained(model).to(device)

    all_embeddings = []
    for i in range(0, len(text), batch_size):
        batch = text[i : i + batch_size]

        inputs = tokenizer(batch, padding=True, truncation=True, return_tensors="pt")
        inputs = {k: v.to(device) for k, v in inputs.items()}

        with torch.no_grad():
            outputs = model(**inputs)

        embeddings = outputs.last_hidden_state[:, 0, :].cpu().tolist()
        all_embeddings.extend(embeddings)

    return all_embeddings


def getDistanceForEmbeddingPair(embeddingA, embeddingB):
    """
    Calculate cosine distance between two lists of embeddings.
    """
    # Normalize the embeddings
    embeddingA_norm = embeddingA / np.linalg.norm(embeddingA, axis=1)[:, np.newaxis]
    embeddingB_norm = embeddingB / np.linalg.norm(embeddingB, axis=1)[:, np.newaxis]

    # Calculate cosine similarity using dot product
    cosine_similarity = np.dot(embeddingB_norm, embeddingA_norm.T)

    # Convert similarity to distance
    cosine_distance = 1 - cosine_similarity

    return cosine_distance


def getClosestMatch(matrix):
    """
    Find closest matching embeddings.
    """

    # Find the index of the minimum distance for each column
    closest_indices = np.argmin(matrix, axis=1)

    # Find the minimum distance values for each column
    min_distances = np.min(matrix, axis=1)

    return np.array([closest_indices, min_distances]).T


def generateSummary(target, disease="inflammatory bowel disease"):

    # Load the pre-trained GPT-2 model and tokenizer
    model_name = "gpt2"
    tokenizer = GPT2Tokenizer.from_pretrained(model_name)
    model = GPT2LMHeadModel.from_pretrained(model_name)

    # Define the prompt for summarization
    prompt = "Summarize in one short paragraph what is known about {} as a therapeutic target for {}.".format(
        target, disease
    )

    # Tokenize the input prompt
    input_ids = tokenizer.encode(prompt, return_tensors="pt")

    # Generate the output using the model
    output_ids = model.generate(
        input_ids,
        max_length=100,  # Limit the length of the generated text
        num_beams=5,  # Use beam search for better results
        early_stopping=True,
    )

    # Decode the output tokens to get the text
    output_text = tokenizer.decode(output_ids[0], skip_special_tokens=True)

    print(output_text)


otDataDir = "../data/ot/"
chemblDataDir = "../data/chembl/"
mixDataDir = "../data/mix/"
ukbbDataDir = "../data/ukbb/"
pathwaysDataDir = "../data/pathways/"

# IBD drug names extracted from UKBB prescriptuons using NER
input = ukbbDataDir + "eidToPrescription.tsv"
output1 = ukbbDataDir + "prescription.tsv"
output2 = ukbbDataDir + "drugName.tsv"
if (not os.path.exists(output1)) or (not os.path.exists(output2)):
    data = nerDrug(
        prescription=pd.read_csv(input, sep="\t", names=["eid", "prescription"])
        .prescription.dropna()
        .unique()
        .tolist()
    )  # more than one drugs can be found for the same prescription

    data = pd.DataFrame(data).dropna().drop_duplicates()
    data[["prescription"]].to_csv(
        output1, sep="\t", index=False, header=False
    )  # redundant prescription for different drugs
    data[["drug"]].to_csv(output2, sep="\t", index=False, header=False)

# IBD UKBB drug names embeddings
input = ukbbDataDir + "drugName.tsv"
output = ukbbDataDir + "drugNameEmb.npy"
if not os.path.exists(output):
    np.save(
        output,
        getEmbeddingForText(
            text=[line.strip("\n").split("\t")[0] for line in open(input, "r")],
        ),
    )

# IBD ChEMBL drug names embeddings
input = chemblDataDir + "drugName.tsv"
output = chemblDataDir + "drugNameEmb.npy"
if not os.path.exists(output):
    np.save(
        output,
        getEmbeddingForText(
            text=[line.strip("\n").split("\t")[0] for line in open(input, "r")],
        ),
    )

# cosine distances between UKBB and ChEMBL IBD drug names embeddings
input1 = chemblDataDir + "drugNameEmb.npy"
input2 = ukbbDataDir + "drugNameEmb.npy"
output = mixDataDir + "chemblDrugNameEmbVsUkbbDrugNameEmb.npy"
if not os.path.exists(output):
    np.save(
        output,
        getDistanceForEmbeddingPair(
            embeddingA=np.load(input1), embeddingB=np.load(input2)
        ),
    )

# closest match in ChEMBL IBD drug names for UKBB IBD drug names
input = mixDataDir + "chemblDrugNameEmbVsUkbbDrugNameEmb.npy"
output = mixDataDir + "closestChemblDrugNameEmbForUkbbDrugNameEmb.npy"
if not os.path.exists(output):
    np.save(output, getClosestMatch(matrix=np.load(input)))

# summary for target + IBD relation
input1 = otDataDir + "novelTargetDetail-EFO_0003767.tsv"
input2 = otDataDir + "target.tsv"
output = mixDataDir + "summaries.txt"
if not os.path.exists(output):
    novel = pd.read_csv(input1, sep="\t")
    novel = novel[novel.novelty >= 0.1].targetSymbol.unique()
    for targetSymbol in novel:
        summary = generateSummary(target=targetSymbol)
        print(summary)
