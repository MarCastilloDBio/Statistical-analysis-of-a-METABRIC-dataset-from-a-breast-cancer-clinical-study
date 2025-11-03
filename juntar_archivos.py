import pandas as pd
import numpy as np

expr = pd.read_csv(r"C:\aEntrega2_programacionbio\brca_metabric\data_mrna_illumina_microarray.txt", sep="\t")
zs   = pd.read_csv(r"C:\aEntrega2_programacionbio\brca_metabric\data_mrna_illumina_microarray_zscores_ref_diploid_samples.txt", sep="\t")
meta = pd.read_csv(r"C:\aEntrega2_programacionbio\brca_metabric\data_clinical_sample.txt", sep="\t", comment="#")


# Tomar 200 muestras al azar de las columnas (excepto las dos primeras: Hugo_Symbol y Entrez_Gene_Id)
np.random.seed(42)
sample_names = np.random.choice(expr.columns[2:], size=200, replace=False)

expr_sub = expr[["Hugo_Symbol", "Entrez_Gene_Id"] + sample_names.tolist()]
zs_sub   = zs[["Hugo_Symbol", "Entrez_Gene_Id"] + sample_names.tolist()]

# Convertir matrices (ancho → largo)
expr_long = expr_sub.melt(
    id_vars=["Hugo_Symbol", "Entrez_Gene_Id"],
    var_name="Sample",
    value_name="Expression"
)

zs_long = zs_sub.melt(
    id_vars=["Hugo_Symbol", "Entrez_Gene_Id"],
    var_name="Sample",
    value_name="Zscore"
)

# Unir expresión + z-scores
df = expr_long.merge(
    zs_long,
    on=["Hugo_Symbol", "Entrez_Gene_Id", "Sample"]
)

# Unir con metadata
df = df.merge(
    meta,
    left_on="Sample",
    right_on="SAMPLE_ID",
    how="inner"  # solo muestras que existen en ambos
)

# Guardar CSV reducido
df.to_csv("METABRIC_brcadata_subset200.csv", index=False)

# Verificar
unique_samples = df["Sample"].unique()
print("Número de muestras únicas en el subset:", len(unique_samples))
