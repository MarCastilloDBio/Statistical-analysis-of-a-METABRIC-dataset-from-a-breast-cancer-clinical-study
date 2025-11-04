import pandas as pd
import random
import seaborn as sns
import matplotlib.pyplot as plt

def GC_content(sequence):
    count_A = sequence.count('A')
    count_G = sequence.count('G')
    count_C = sequence.count('C')
    count_T = sequence.count('T')
    GC_content = (count_G + count_C) / len(sequence)
    GC_content = round(100 * GC_content, 2)
    return GC_content
    pass

filename ='ecoli_seq.fasta'
sequence = str()
with open(filename) as file:
    next(file)
    for line in file:
        line = line.strip()
        sequence += line

sequence = list(sequence)
random.shuffle(sequence)
sequence = "".join(sequence)[:100000]

gc_list = []
seqlen = 100
for repetition in range(1000):
    random_seq = list(sequence)
    random.shuffle(random_seq)
    random_seq = "".join(random_seq)[:seqlen]
    gc_value = GC_content(random_seq)
    gc_list.append(gc_value)

gc_value = pd.Series(gc_list)

# Calcular valores p para ambos extremos
p_value_superior = sum(gc_value >= 60) / repetition
p_value_inferior = sum(gc_value <= 40) / repetition  # umbral inferior

# 4. Gráfico 1: Histograma de la simulación con ambos Valores P
plt.figure(figsize=(10, 6))
sns.histplot(gc_list, kde=True, bins=30, color='skyblue', edgecolor='black')

# Marcar el umbral superior
plt.axvline(60, color='red', linestyle='--', linewidth=2, label=f'Umbral superior GC% ({60}%)')

# Marcar el umbral inferior
plt.axvline(40, color='orange', linestyle='--', linewidth=2, label=f'Umbral inferior GC% ({40}%)')

# Añadir el valor p superior
plt.text(0.95, 0.95,
         f'$P_{{empirical}}$ ($\geq$ {60}% GC) = {p_value_superior:.4f}',
         transform=plt.gca().transAxes,
         fontsize=12,
         verticalalignment='top',
         horizontalalignment='right',
         bbox=dict(boxstyle="round,pad=0.5", fc="white", alpha=0.7))

# Añadir el valor p inferior
plt.text(0.95, 0.85,
         f'$P_{{empirical}}$ ($\leq$ {40}% GC) = {p_value_inferior:.4f}',
         transform=plt.gca().transAxes,
         fontsize=12,
         verticalalignment='top',
         horizontalalignment='right',
         bbox=dict(boxstyle="round,pad=0.5", fc="white", alpha=0.7))

plt.title(f'Distribución de Contenido GC')
plt.xlabel('Contenido GC (%)')
plt.ylabel('Frecuencia')
plt.legend()
plt.grid(axis='y', alpha=0.5)
plt.show()

