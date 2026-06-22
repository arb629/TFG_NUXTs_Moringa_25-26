import pandas as pd

cols = ["qseqid","sseqid","pident","length","mismatch","gapopen",
        "qstart","qend","sstart","send","evalue","bitscore"]

def get_best_hits(blast_file):
    df = pd.read_csv(blast_file, sep="\t", header=None, names=cols)
    return df.sort_values("bitscore", ascending=False).drop_duplicates("qseqid").set_index("qseqid")["sseqid"]

for tipo in ["nupt", "numt"]:
    fwd = get_best_hits(f"/home/kubuntu/Andrea/TFG/Diagrama_venn/{tipo}_oleif_vs_steno.out")
    rev = get_best_hits(f"/home/kubuntu/Andrea/TFG/Diagrama_venn/{tipo}_steno_vs_oleif.out")

    # RBH: oleif->steno y steno->oleif son mejores hits mutuos
    rbh = set()
    for oleif_id, steno_id in fwd.items():
        if steno_id in rev.index and rev[steno_id] == oleif_id:
            rbh.add(oleif_id)

    all_oleif = set(fwd.index)
    all_steno = set(rev.index)
    comunes   = rbh
    # Guardar también los IDs de stenopetala correspondientes a los comunes
    comunes_steno = {steno_id for steno_id, oleif_id in rev.items() if oleif_id in comunes}
    esp_oleif = all_oleif - comunes
    esp_steno = all_steno - {rev[s] for s in all_steno if rev[s] in comunes}

    print(f"\n{tipo.upper()}:")
    print(f"  Total oleif:       {len(all_oleif)}")
    print(f"  Total steno:       {len(all_steno)}")
    print(f"  Comunes (RBH):     {len(comunes)}")
    print(f"  Especificos oleif: {len(esp_oleif)}")
    print(f"  Especificos steno: {len(all_steno - set(rev[rev.isin(comunes)].index))}")

    # Guardar listas
    pd.Series(list(comunes)).to_csv(f"/home/kubuntu/Andrea/TFG/Diagrama_venn/{tipo}_comunes.txt", index=False, header=False)
    pd.Series(list(esp_steno)).to_csv(f"/home/kubuntu/Andrea/TFG/Diagrama_venn/{tipo}_especificos_steno.txt", index=False, header=False)
    pd.Series(list(esp_oleif)).to_csv(f"/home/kubuntu/Andrea/TFG/Diagrama_venn/{tipo}_especificos_oleif.txt", index=False, header=False)
    pd.Series(list(comunes_steno)).to_csv(f"/home/kubuntu/Andrea/TFG/Diagrama_venn/{tipo}_comunes_steno.txt", index=False, header=False)
