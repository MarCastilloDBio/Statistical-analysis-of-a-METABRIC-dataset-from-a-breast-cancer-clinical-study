import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import norm
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

# --------------------------
# Parámetros
CSV_FILE = r"C:\aEntrega2_programacionbio\METABRIC_brcadata_subset200.csv"
SAMPLE_ID = "MB-4738"
Z_CUTOFF = 2.0         # ±2 para decisión
NUM_LABELS = 12        # cuantos genes extremos anotar

sns.set(style="whitegrid")

# --------------------------
# 1) Cargar y filtrar
df = pd.read_csv(CSV_FILE)
df_s = df[df["Sample"] == SAMPLE_ID].copy()

# --------------------------
# 2) Z-score a p (two-tailed) y eje Y
df_s["p_value"] = norm.sf(np.abs(df_s["Zscore"])) * 2.0
# Evitar -log10(0)
df_s["p_value"] = df_s["p_value"].clip(lower=1e-300)
df_s["neg_log10_p"] = -np.log10(df_s["p_value"])

# --------------------------
# 3) Preparar mapa cromatico
# Rango de z a considerar para escala (más allá de +/-2)
max_pos = max(df_s["Zscore"].max(), Z_CUTOFF + 0.1)
max_neg = min(df_s["Zscore"].min(), -Z_CUTOFF - 0.1)

# Normalizadores: mapean Z a [0,1] para el cmap
norm_pos = Normalize(vmin=Z_CUTOFF, vmax=max_pos)
norm_neg = Normalize(vmin=Z_CUTOFF, vmax=abs(max_neg))  # se usa abs para normalizar magnitud

cmap_red = plt.cm.Reds       # claro->oscuro cuando el valor del norm aumenta
cmap_blue = plt.cm.Blues     # claro->oscuro

# Asignar color por fila
colors = []
for z in df_s["Zscore"]:
    if z >= Z_CUTOFF:
        # valor normalizado desde Z_CUTOFF (claro) hasta max_pos (oscuro)
        v = norm_pos(z)
        colors.append(cmap_red(v))
    elif z <= -Z_CUTOFF:
        v = norm_neg(abs(z))
        colors.append(cmap_blue(v))
    else:
        colors.append((0.7, 0.7, 0.7, 0.9))  # gris claro, alpha 0.9

df_s["color"] = colors

# --------------------------
# 4) Gráfico con seaborn/matplotlib
fig, ax = plt.subplots(figsize=(14, 9))

# primero se trazan todos los puntos (usa ax.scatter para colores ya calculados)
ax.scatter(df_s["Zscore"], df_s["neg_log10_p"], c=df_s["color"], s=28, edgecolor='k', linewidth=0.2, alpha=0.85)

# Líneas de corte
ax.axvline(x= Z_CUTOFF, color='black', linestyle='--', linewidth=1)
ax.axvline(x=-Z_CUTOFF, color='black', linestyle='--', linewidth=1)
# línea horizontal para p = 0.05
ax.axhline(y=-np.log10(0.05), color='black', linestyle='--', linewidth=1)

# Etiquetas y formato
ax.set_title(f"Scatter plot — Muestra {SAMPLE_ID}", fontsize=16, weight='bold')
ax.set_xlabel("Log2-transformed mRNA relative z-Scores", fontsize=13)
ax.set_ylabel("-log10(p-value) (two-tailed)", fontsize=13)

# Grid suave
ax.grid(True, linestyle=':', alpha=0.4)

# --------------------------
# 5) Colorbars (una para rojo, otra para azul)
# Red colorbar (sobreexpresión): map from Z_CUTOFF to max position
sm_red = ScalarMappable(norm=norm_pos, cmap=cmap_red)
sm_red.set_array([])

# Blue colorbar (subexpresión): map from -Z_CUTOFF yo min position (se usa magnitud)
sm_blue = ScalarMappable(norm=norm_neg, cmap=cmap_blue)
sm_blue.set_array([])

# Añadir colorbars a la figura
cbar_red = fig.colorbar(sm_red, ax=ax, fraction=0.046, pad=0.06)
cbar_red.set_label('Sobreexpresión (Z-score)', fontsize=12)
cbar_blue = fig.colorbar(sm_blue, ax=ax, fraction=0.046, pad=0.11)
cbar_blue.set_label('Subexpresión (abs Z-score)', fontsize=12)

# --------------------------
# 6) Anotar genes más extremos
df_s["absZ"] = df_s["Zscore"].abs()
top = df_s.nlargest(NUM_LABELS, "absZ")

# Para evitar superposición simple: alternar posiciones arriba/abajo y aplicar pequeño desplazamiento x
y_offset = (df_s["neg_log10_p"].max() - df_s["neg_log10_p"].min()) * 0.02
for i, (_, row) in enumerate(top.iterrows()):
    x = row["Zscore"]
    y = row["neg_log10_p"]
    # alternar etiqueta hacia arriba o abajo
    dy = y_offset if (i % 2 == 0) else -y_offset*1.5
    dx = 0.02 * x if x != 0 else 0.1
    ax.annotate(row["Hugo_Symbol"], xy=(x, y), xytext=(x+dx, y+dy),
                fontsize=9, weight='bold',
                arrowprops=dict(arrowstyle='-', linewidth=0.5, color='gray', alpha=0.6))

# --------------------------
# 7) Leyendas
ax.scatter([], [], c='lightgrey', s=30, edgecolor='k', label=f'No significativo (|Z| < {Z_CUTOFF})')
ax.scatter([], [], c='red', s=30, edgecolor='k', label=f'Sobreexpresión (Z >= {Z_CUTOFF})')
ax.scatter([], [], c='blue', s=30, edgecolor='k', label=f'Subexpresión (Z <= -{Z_CUTOFF})')
ax.legend(loc='upper right', frameon=True)

plt.tight_layout()
plt.show()


# r"C:\aEntrega2_programacionbio\METABRIC_brcadata_subset200.csv"