import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

# ------------------------------
# Cargar datos
df = pd.read_csv(r"C:\aEntrega2_programacionbio\METABRIC_brcadata_subset200.csv")
# ------------------------------
# Seleccionar genes
significant_genes = [
    "ASXL3", "C20orf144", "ALKAL2", "TNNNSL6", "TMEM178A",
    "LAMC3", "SMYD1", "OXTR", "ARHGAP36", "NEUROG2"
]
df_sig = df[df['Hugo_Symbol'].isin(significant_genes)].copy()
# ------------------------------
# Crear matriz genes x muestras
heatmap_data = df_sig.pivot(index="Hugo_Symbol", columns="Sample", values="Zscore")
heatmap_data = heatmap_data.fillna(0)  # reemplazar NA por 0
# ------------------------------
# Crear barras de color para cada subtipo de cancer
sample_types = df_sig[['Sample', 'CANCER_TYPE_DETAILED']].drop_duplicates().set_index('Sample')
subtype_colors = {
    "Breast Mixed Ductal and Lobular Carcinoma": "#1f77b4",  # azul oscuro
    "Breast Invasive Ductal Carcinoma": "#ff7f0e",           # naranja
    "Breast Invasive Lobular Carcinoma": "#2ca02c",          # verde
}
col_colors = sample_types['CANCER_TYPE_DETAILED'].map(subtype_colors)
# ------------------------------
# Crear gradiente
colors = ["#FFA500", "#FFB84D", "#FFCC99", "#E6CCE6", "#CC99CC", "#9966CC", "#800080"]
cmap = LinearSegmentedColormap.from_list("orange_to_purple", colors, N=256)
# ------------------------------
# Configurar tema
sns.set_theme(style="white")
# ------------------------------
# Crear clustermap
cg = sns.clustermap(
    heatmap_data,
    cmap=cmap,  # gradiente naranja → morado oscuro
    standard_scale=0,
    linewidths=0.5,
    figsize=(16,10),
    col_cluster=True,
    row_cluster=True,
    col_colors=col_colors
)
# ------------------------------
# Cambiar etiquetas de los ejes
cg.ax_heatmap.set_ylabel('Genes', fontsize=12)
cg.ax_heatmap.set_xlabel('Muestras', fontsize=12)
# ------------------------------
# Leyenda lateral de los subtipos
for label in subtype_colors:
    cg.ax_col_dendrogram.bar(0, 0, color=subtype_colors[label], label=label, linewidth=0)
cg.ax_col_dendrogram.legend(title="Subtipo de cáncer", loc="center", ncol=1)
plt.suptitle("Clustermap — Genes alterados (METABRIC)", fontsize=16)
plt.show()


