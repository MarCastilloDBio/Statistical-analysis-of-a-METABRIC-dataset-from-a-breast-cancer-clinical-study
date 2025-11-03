import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# CARGAR DATOS
df = pd.read_csv(r"C:\aEntrega2_programacionbio\METABRIC_brcadata_subset200.csv")

# GENES CON MAYOR SUB Y SOBREXPRESIÓN
significant_genes = [
    "ASXL3", "C20orf144", "ALKAL2", "TNNNSL6", "TMEM178A",
    "LAMC3", "SMYD1", "OXTR", "ARHGAP36", "NEUROG2"
]

# FILTRAR SOLO ESOS GENES
df_sig = df[df["Hugo_Symbol"].isin(significant_genes)]

# CREAR MATRIZ GENES x MUESTRAS
heatmap_data = df_sig.pivot_table(
    index="Hugo_Symbol",
    columns="Sample",
    values="Zscore"
)

# ORDENAR GENES
heatmap_data = heatmap_data.loc[
    heatmap_data.mean(axis=1).sort_values(ascending=False).index
]

# COLORMAP MÁS INTENSO - se usan tonos de colores más saturados
# Azul oscuro → Azul → Blanco → Rojo → Rojo oscuro
strong_cmap = LinearSegmentedColormap.from_list(
    "strongmap",
    ["#0011FF", "#0055FF", "white", "#FF4444", "#B30000"]
)

# HEATMAP
plt.figure(figsize=(22, 6))

sns.heatmap(
    heatmap_data,
    cmap=strong_cmap,  # mapa más intenso
    center=0,
    linewidths=0.2,
    linecolor="gray",
    cbar_kws={"label": "Log2 Z-score"}
)

plt.title("Heatmap — Genes significativamente alterados (METABRIC Breast Cancer)")
plt.xlabel("Muestras")
plt.ylabel("Genes")
plt.xticks(rotation=90)
plt.yticks(rotation=0)
plt.tight_layout()
plt.show()