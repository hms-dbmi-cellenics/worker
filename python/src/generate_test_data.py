# This script generates test data for working on data sets without the backend.

import anndata
import scanpy as sc
import json

adata = anndata.read("./basic-filtered.h5ad")

# Do UMAP and louvain
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.tl.louvain(adata)

print(adata)
sc.tl.dendrogram(adata, groupby="louvain")
print(adata.uns["dendrogram_louvain"])


# Get a list of cells & cluster info
cells = adata.obs.index.tolist()
clusters = adata.obs["louvain"].tolist()
embedding = adata.obsm["X_umap"].tolist()

# Categorical embedding schema
embedding_categorical = {
    "cells": cells,
    "categories": [
        {"categoryName": "louvain", "type": "categorical", "values": clusters}
    ],
    "embedding": embedding,
}

# Get expression level of an interesting gene (CST3)
interesting_gene = adata.raw[:, adata.raw.var.index == "CST3"]
expression_level = interesting_gene.X.toarray().flatten().tolist()
embedding_continuous = {
    "cells": cells,
    "categories": [
        {"categoryName": "CST3", "type": "continuous", "values": expression_level}
    ],
    "embedding": embedding,
}

# Get heatmap data
adata = anndata.read("./basic-filtered.h5ad")
adata = adata.raw.to_adata()
adata = adata[:, adata.var.highly_variable]

gene_indices = adata.var.sort_values(by="dispersions", ascending=False).index
gene_indices = gene_indices[0:250].tolist()


heatmap = {
    "cellNames": cells,
    "categories": [
        {"categoryName": "louvain", "type": "categorical", "values": clusters}
    ],
    "heatmapData": [],
}


for gene in gene_indices:
    view = adata[:, adata.var.index == gene]
    expression = view.X.toarray().flatten().tolist()

    heatmap["heatmapData"].append({"geneName": gene, "expression": expression})

with open("heatmap.json", "w") as f:
    json.dump(heatmap, f)

with open("embedding_categorical.json", "w") as f:
    json.dump(embedding_categorical, f)

with open("embedding_continuous.json", "w") as f:
    json.dump(embedding_continuous, f)
