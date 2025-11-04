import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# ------------------------------
# Cargar datos
# ------------------------------
df = pd.read_csv(r"C:\aEntrega2_programacionbio\METABRIC_brcadata_subset200.csv")

# ------------------------------
# Definir categorías combinando ER, PR, HER2
# ------------------------------
conditions = [
    (df['ER_STATUS'] == 'Positive') & (df['PR_STATUS'] == 'Positive') & (df['HER2_STATUS'] == 'Negative'),
    (df['ER_STATUS'] == 'Negative') & (df['PR_STATUS'] == 'Negative') & (df['HER2_STATUS'] == 'Negative'),
    (df['HER2_STATUS'] == 'Positive')
]
choices = ['Luminal-like', 'Triple-Negative', 'HER2-Enriched']
df['Subtype'] = np.select(conditions, choices, default='Other')

# ------------------------------
# Filtrar solo las categorías de interés
# ------------------------------
df_plot = df[df['Subtype'] != 'Other'].copy()

# ------------------------------
# Configurar tema de Seaborn
# ------------------------------
sns.set_theme(style="whitegrid")
palette_colors = {"Luminal-like": "green", "Triple-Negative": "blue", "HER2-Enriched": "violet"}

# ------------------------------
# Crear JointGrid para más control
# ------------------------------
g = sns.JointGrid(data=df_plot, x="TUMOR_SIZE", y="TMB_NONSYNONYMOUS", height=8, space=0.2)

# Scatter central por subtipo
for subtype, color in palette_colors.items():
    subset = df_plot[df_plot['Subtype'] == subtype]
    g.ax_joint.scatter(subset['TUMOR_SIZE'], subset['TMB_NONSYNONYMOUS'],
                       c=color, label=subtype, alpha=0.6, s=50)

# Histogramas marginales con KDE por subtipo
for subtype, color in palette_colors.items():
    subset = df_plot[df_plot['Subtype'] == subtype]
    sns.kdeplot(data=subset, x='TUMOR_SIZE', ax=g.ax_marg_x, color=color, fill=True, alpha=0.3)
    sns.kdeplot(data=subset, y='TMB_NONSYNONYMOUS', ax=g.ax_marg_y, color=color, fill=True, alpha=0.3)

# Leyenda y etiquetas
g.ax_joint.legend(title="Subtipo", loc="upper right")
g.set_axis_labels("Tamaño del tumor (mm)", "Carga mutacional tumoral (TMB)", fontsize=12)

# Título con más espacio
g.fig.suptitle("Relación entre tamaño del tumor y TMB según subtipo (METABRIC)",
               fontsize=14, y=0.98)

plt.tight_layout()
plt.show()







