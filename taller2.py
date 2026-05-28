import pandas as pd
import numpy as np

# --- Carga de Datos ---
# Cargamos el CSV 
df = pd.read_csv('datos_genotipo.csv', index_col=0)
G = df.values # Matriz G
n, m = G.shape # n=30 individuos, m=100 SNPs

# --- Estimación de Frecuencias y Limpieza ---
# Calculamos la frecuencia p_j para cada SNP (columna)
p_j = np.sum(G, axis=0) / (2 * n)

# Filtramos singularidades
# Conservamos solo las columnas donde la frecuencia no sea 0 ni 1
valid_snps_mask = (p_j > 0) & (p_j < 1)

# Aplicamos la máscara para obtener la matriz y frecuencias limpias
G_clean = G[:, valid_snps_mask]
p_j_clean = p_j[valid_snps_mask]

# Mostramos los resultados en consola para validar
print(f"Dimensiones originales de G: {G.shape}")
print(f"Dimensiones de G despues de la limpieza: {G_clean.shape}")
print(f"SNPs (columnas) descartados por ser constantes: {m - G_clean.shape[1]}")