import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Ellipse

mpl.rcParams["text.color"] = "black"
mpl.rcParams["axes.labelcolor"] = "black"
mpl.rcParams["xtick.color"] = "black"
mpl.rcParams["ytick.color"] = "black"

datos = {
    "NUPTs": {"oleif": 2814, "comunes": 196, "steno": 776},
    "NUMTs": {"oleif": 3130, "comunes": 186, "steno": 336}
}

colores = {
    "NUPTs": {"oleif": "#4f7a4d", "steno": "#7baf48"},
    "NUMTs": {"oleif": "#a13f2b", "steno": "#df3f2d"}
}

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.patch.set_facecolor("white")

for ax, (tipo, vals) in zip(axes, datos.items()):
    col_oleif = colores[tipo]["oleif"]
    col_steno = colores[tipo]["steno"]

    ax.set_xlim(0, 10)
    ax.set_ylim(0, 7)
    ax.set_aspect("equal")
    ax.axis("off")

    ell_oleif = Ellipse(
        xy=(3.8, 3.5), width=5.8, height=4.5,
        facecolor=col_oleif, alpha=0.4,
        edgecolor=col_oleif, linewidth=2.5,
        zorder=2, antialiased=True
    )
    ell_steno = Ellipse(
        xy=(6.2, 3.5), width=5.8, height=4.5,
        facecolor=col_steno, alpha=0.3,
        edgecolor=col_steno, linewidth=2.5,
        zorder=2, antialiased=True
    )

    ax.add_patch(ell_oleif)
    ax.add_patch(ell_steno)

    from matplotlib.patches import Ellipse
    borde_oleif = Ellipse(
        xy=(3.8, 3.5), width=5.8, height=4.5,
        facecolor="none",
        edgecolor=col_oleif, linewidth=2.5, zorder=3
    )
    borde_steno = Ellipse(
        xy=(6.2, 3.5), width=5.8, height=4.5,
        facecolor="none",
        edgecolor=col_steno, linewidth=2.5, zorder=3
    )
    ax.add_patch(borde_oleif)
    ax.add_patch(borde_steno)

    ax.text(2.2, 3.5, str(vals["oleif"]),
            ha="center", va="center", fontsize=18,
            fontweight="bold", color=col_oleif)

    ax.text(5.0, 3.5, str(vals["comunes"]),
            ha="center", va="center", fontsize=16,
            fontweight="bold", color="#444444")

    ax.text(7.8, 3.5, str(vals["steno"]),
            ha="center", va="center", fontsize=18,
            fontweight="bold", color=col_steno)

    ax.text(2.8, 1.0, r"$\it{M. oleifera}$",
            ha="center", va="center", fontsize=13,
            fontweight="bold", color=col_oleif, zorder=5)

    ax.text(7.2, 1.0, r"$\it{M. stenopetala}$",
            ha="center", va="center", fontsize=13,
            fontweight="bold", color=col_steno, zorder=5)

    n_total_oleif = vals["oleif"] + vals["comunes"]
    n_total_steno = vals["steno"] + vals["comunes"]

    ax.text(2.8, 0.45, f"n = {n_total_oleif}",
            ha="center", va="center", fontsize=11, color="black")

    ax.text(7.2, 0.45, f"n = {n_total_steno}",
            ha="center", va="center", fontsize=11, color="black")

    ax.set_title(tipo, fontsize=16, fontweight="bold", color="black", pad=10)

fig.suptitle("NUXTs específicos y comunes entre especies",
             fontsize=14, fontweight="bold", color="#1a3a5c", y=1.02)

plt.tight_layout()
fig.savefig("venn_NUXTs_v3.png", dpi=600, bbox_inches="tight",
            facecolor=fig.get_facecolor(), edgecolor="none")
plt.show()
