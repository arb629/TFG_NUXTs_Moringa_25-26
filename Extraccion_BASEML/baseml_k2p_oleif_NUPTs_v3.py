# ============================================================
# BASEML K2P — M. oleifera NUPTs filtrados (sin IR)  v3
# Input:  oleif_ptld_nuc_CORREGIDO.out
# Output: /home/kubuntu/Andrea/TFG/Filtro_Blast_IR/v3/oleif/NUPTs/
# Cambios v3:
#   - Filtro NOT CONVERGENT
#   - outfile en .ctl con ruta absoluta
# ============================================================

from Bio import SeqIO
import subprocess
import os

# ── RUTAS ────────────────────────────────────────────────────
out_file       = "/home/kubuntu/Andrea/TFG/M_oleifera/NUXTs/oleif_ptld_nuc_CORREGIDO.out"
nuclear_fasta  = "/home/kubuntu/Andrea/TFG/M_oleifera/NU/MoringaV2.genome.fa"
organelo_fasta = "/home/kubuntu/Andrea/TFG/M_oleifera/PT/m_oleifera_PT.fasta"
output_dir     = "/home/kubuntu/Andrea/TFG/Filtro_Blast_IR/v3/oleif/NUPTs"
os.makedirs(output_dir, exist_ok=True)

# ── CARGAR GENOMAS ───────────────────────────────────────────
print("Cargando genomas...")
nuclear      = {rec.id: rec for rec in SeqIO.parse(nuclear_fasta,  "fasta")}
organelo_raw = {rec.id: rec for rec in SeqIO.parse(organelo_fasta, "fasta")}
mapeo_organelo = {"ctg000001c": "m_oleifera_PT"}
organelo = {}
for nombre_out, nombre_fasta in mapeo_organelo.items():
    if nombre_fasta in organelo_raw:
        organelo[nombre_out] = organelo_raw[nombre_fasta]
print(f"Nuclear:  {len(nuclear)} secuencias")
print(f"Organelo: {len(organelo)} secuencias")

procesados    = 0
descartados   = 0
no_convergent = 0


def check_convergence(baseml_out_path):
    if not os.path.exists(baseml_out_path) or os.path.getsize(baseml_out_path) == 0:
        return False
    with open(baseml_out_path, "r") as f:
        content = f.read()
    if "NOT CONVERGENT" in content:
        return False
    return True


with open(out_file) as f:
    for i, line in enumerate(f):
        line = line.strip()
        if not line:
            continue
        cols = line.split("\t")
        if len(cols) < 12:
            continue

        org_name  = cols[0]
        nuc_chr   = cols[1]
        org_start = int(cols[6])
        org_end   = int(cols[7])
        nuc_start = int(cols[8])
        nuc_end   = int(cols[9])

        org_s = min(org_start, org_end) - 1
        org_e = max(org_start, org_end)
        nuc_s = min(nuc_start, nuc_end) - 1
        nuc_e = max(nuc_start, nuc_end)

        if org_name not in organelo or nuc_chr not in nuclear:
            continue

        seq_org = organelo[org_name].seq[org_s:org_e]
        seq_nuc = nuclear[nuc_chr].seq[nuc_s:nuc_e]

        if len(seq_org) < 100 or len(seq_nuc) < 100:
            continue

        nombre     = f"par_{i+1}_pt"
        fasta_file = os.path.join(output_dir, f"{nombre}.fasta")
        aln_file   = os.path.join(output_dir, f"{nombre}_aln.fasta")
        trim_file  = os.path.join(output_dir, f"{nombre}_trim.fasta")
        phy_file   = os.path.join(output_dir, f"{nombre}.phy")
        tree_file  = os.path.join(output_dir, f"{nombre}.tree")
        ctl_k2p    = os.path.join(output_dir, f"{nombre}_K2P.ctl")
        out_k2p    = os.path.join(output_dir, f"{nombre}_K2P.out")

        with open(fasta_file, "w") as out:
            out.write(f">nuclear_{nuc_chr}_{nuc_s}_{nuc_e}\n{seq_nuc}\n")
            out.write(f">organelo_{org_s}_{org_e}\n{seq_org}\n")

        subprocess.run(
            f"mafft --auto {fasta_file} > {aln_file} 2>/dev/null",
            shell=True)

        subprocess.run(
            f"trimal -in {aln_file} -out {trim_file} -automated1 -fasta",
            shell=True, capture_output=True, text=True)

        if not os.path.exists(trim_file) or os.path.getsize(trim_file) == 0:
            descartados += 1
            continue

        records_trim = list(SeqIO.parse(trim_file, "fasta"))
        if len(records_trim) < 2:
            descartados += 1
            continue

        longitud_trim = len(records_trim[0].seq)
        if longitud_trim < 50:
            descartados += 1
            continue

        with open(phy_file, "w") as out:
            out.write(f" 2 {longitud_trim}\n")
            for rec in records_trim:
                name = rec.id[:10].ljust(10)
                out.write(f"{name}  {str(rec.seq)}\n")

        with open(tree_file, "w") as out:
            out.write("1\n(nuclear, organelo);\n")

        with open(ctl_k2p, "w") as out:
            out.write(f"seqfile  = {phy_file}\n")
            out.write(f"treefile = {tree_file}\n")
            out.write(f"outfile  = {out_k2p}\n")
            out.write("noisy = 0\nverbose = 0\nrunmode = 0\n")
            out.write("model = 1\nMgene = 0\nfix_kappa = 0\nkappa = 2\n")
            out.write("fix_alpha = 1\nalpha = 0\ngetSE = 0\nRateAncestor = 0\n")
            out.write("Small_Diff = 1e-8\n")
            out.write("cleandata = 1\n")

        subprocess.run(
            f"cd {output_dir} && baseml {ctl_k2p}",
            shell=True, capture_output=True)

        # ── FILTRO NOT CONVERGENT ────────────────────────────
        if not check_convergence(out_k2p):
            no_convergent += 1
            descartados += 1
            continue

        procesados += 1
        if procesados % 100 == 0:
            print(f"Procesados: {procesados} | Descartados: {descartados} "
                  f"(no convergent: {no_convergent})")

print(f"\nDONE")
print(f"  Procesados (válidos):    {procesados}")
print(f"  Descartados (total):     {descartados}")
print(f"  → Por no convergencia:   {no_convergent}")
