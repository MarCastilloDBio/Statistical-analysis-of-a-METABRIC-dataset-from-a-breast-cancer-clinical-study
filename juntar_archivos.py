import pandas as pd
import numpy as np

expr = pd.read_csv(r"C:\aEntrega2_programacionbio\brca_metabric\data_mrna_illumina_microarray.txt", sep="\t")
zs   = pd.read_csv(r"C:\aEntrega2_programacionbio\brca_metabric\data_mrna_illumina_microarray_zscores_ref_diploid_samples.txt", sep="\t")
meta = pd.read_csv(r"C:\aEntrega2_programacionbio\brca_metabric\data_clinical_sample.txt", sep="\t", comment="#")


# 2️⃣ Convertir matrices (ancho → largo)
expr_long = expr.melt(
    id_vars=["Hugo_Symbol", "Entrez_Gene_Id"],
    var_name="Sample",
    value_name="Expression"
)

zs_long = zs.melt(
    id_vars=["Hugo_Symbol", "Entrez_Gene_Id"],
    var_name="Sample",
    value_name="Zscore"
)

# 3️⃣ Unir expresión + z-scores
df = expr_long.merge(
    zs_long,
    on=["Hugo_Symbol", "Entrez_Gene_Id", "Sample"]
)

# 4️⃣ Unir con metadata
df = df.merge(
    meta,
    left_on="Sample",
    right_on="SAMPLE_ID",
    how="inner"  # solo muestras que existen en ambos
)

# Escoge 1000 SAMPLE_ID únicos al azar
np.random.seed(42)  # para reproducibilidad
sample_subset = np.random.choice(df['Sample'].unique(), size=500, replace=False)

# Filtrar el dataframe
df_subset = df[df['Sample'].isin(sample_subset)]

# 6️⃣ Guardar CSV reducido
df_subset.to_csv("METABRIC_brcadata_subset500.csv", index=False)

print("Archivo CSV con 1000 muestras generado con éxito ✅")