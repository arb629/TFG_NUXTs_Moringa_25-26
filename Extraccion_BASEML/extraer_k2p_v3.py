# ============================================================
# Extraer valores K2P de BASEML — v3
# Lee los _K2P.out de las 4 carpetas y genera 4 CSVs
# Output: /home/kubuntu/Andrea/TFG/Filtro_Blast_IR/v3/
# ============================================================

import glob
import re
import os

combinaciones = [
    (
        "oleif NUPTs",
        "/home/kubuntu/Andrea/TFG/Filtro_Blast_IR/v3/oleif/NUPTs",
        "/home/kubuntu/Andrea/TFG/Filtro_Blast_IR/v3/k_NUPTs_oleif_v3.csv",
    ),
    (
        "oleif NUMTs",
        "/home/kubuntu/Andrea/TFG/Filtro_Blast_IR/v3/oleif/NUMTs",
        "/home/kubuntu/Andrea/TFG/Filtro_Blast_IR/v3/k_NUMTs_oleif_v3.csv",
    ),
    (
        "steno NUPTs",
        "/home/kubuntu/Andrea/TFG/Filtro_Blast_IR/v3/steno/NUPTs",
        "/home/kubuntu/Andrea/TFG/Filtro_Blast_IR/v3/k_NUPTs_steno_v3.csv",
    ),
    (
        "steno NUMTs",
        "/home/kubuntu/Andrea/TFG/Filtro_Blast_IR/v3/steno/NUMTs",
        "/home/kubuntu/Andrea/TFG/Filtro_Blast_IR/v3/k_NUMTs_steno_v3.csv",
    ),
]

for nombre, input_dir, output_csv in combinaciones:
    print(f"\nProcesando {nombre}...")
    resultados = []
    archivos = sorted(glob.glob(os.path.join(input_dir, "*_K2P.out")))

    if not archivos:
        print(f"  AVISO: no se encontraron archivos _K2P.out en {input_dir}")
        continue

    for out_file in archivos:
        match_par = re.search(r'par_(\d+)', out_file)
        if not match_par:
            continue
        par = match_par.group(1)

        with open(out_file) as f:
            contenido = f.read()

        # Saltar archivos sin resultado o con no-convergencia (doble filtro)
        if "Distances" not in contenido:
            continue
        if "NOT CONVERGENT" in contenido:
            continue

        match_k = re.search(
            r'\n\S+\s+([\d\.]+)[\(\s]',
            contenido.split("Distances")[1]
        )
        if match_k:
            resultados.append((int(par), "K2P", match_k.group(1)))

    resultados.sort(key=lambda x: x[0])

    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    with open(output_csv, "w") as f:
        f.write("Fragmento\tModelo\tk\n")
        for par, modelo, k in resultados:
            f.write(f"{par}\t{modelo}\t{k}\n")

    print(f"  Pares extraídos: {len(resultados)}")
    print(f"  CSV guardado en: {output_csv}")

print("\nDONE")
