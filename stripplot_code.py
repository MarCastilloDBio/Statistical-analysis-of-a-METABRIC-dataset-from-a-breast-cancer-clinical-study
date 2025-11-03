import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# CARGAR DATOS
df = pd.read_csv(r"C:\aEntrega2_programacionbio\METABRIC_brcadata_subset200.csv")

# FILTRAR SOLO EL GEN ARHGAP36
df_gene = df[df["Hugo_Symbol"] == "ARHGAP36"]

# ESTILO SEABORN
sns.set(style="whitegrid")

plt.figure(figsize=(14, 6))

# BOXPLOT POR SUBTIPO DE CÁNCER DE MAMA
sns.boxplot(
    x="CANCER_TYPE_DETAILED",
    y="Expression",
    data=df_gene,
    whis=1.5,
    width=0.5,
    showcaps=True,
    boxprops={'facecolor':'None'},  #  establecer caja transparente
    showfliers=False,               # no duplicar outliers
    medianprops={'color':'black', 'linewidth':2}
)

# STRIPPLOT SOBRE EL BOXPLOT
sns.stripplot(
    x="CANCER_TYPE_DETAILED",
    y="Expression",
    data=df_gene,
    hue="CANCER_TYPE_DETAILED",   # puntos de diferente color por subtipo
    dodge=False,
    palette="Set2",
    size=6,
    alpha=0.8
)

# AJUSTES DEL GRAFICO
plt.title("Expresión de ARHGAP36 según subtipo de cáncer de mama (METABRIC)")
plt.xlabel("Subtipo de cáncer de mama")
plt.ylabel("Expression log2 intensity levels")
plt.xticks(rotation=45)

plt.tight_layout()
plt.show()