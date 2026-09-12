# Figure 5 — scripts and vector panels

All SVGs in `svg/` are written with `matplotlib.rcParams["svg.fonttype"]="none"`, so every
label is a real `<text>` element and remains an editable text layer in Illustrator/Inkscape.
(Matplotlib writes the font as a stack beginning with DejaVu Sans; if that font is absent
locally the labels will substitute and reflow slightly — reassign the font once on import.)

## Panels

| Figure | script | output |
|---|---|---|
| 5A | *provided by author* — example Hi-C maps + CRUSH tracks | — |
| 5B | `Fig5B_compartment_by_GC.py` | `Fig5B_compartment_by_GC.svg` |
| 5C | *provided by author* — example CTCF loop | — |
| 5D | `Fig5D_APA_all_CTCF_loops.py` | `Fig5D_APA_all_CTCF_loops.svg` |
| 5E | `Fig5E_APA_convergent_GGG.py` | `Fig5E_APA_convergent_GGG.svg` |
| 5F | `Fig5F_anchor_depletion.py` | `Fig5F_anchor_depletion.svg` |
| 5G | `Fig5G_EP_by_GC_MA.py` | `Fig5G_EP_by_GC_MA.svg` |
| 5H | `Fig5H_EP_by_CTCFclass.py` | `Fig5H_EP_by_CTCFclass_MA.svg` **and** `Fig5H_EP_by_CTCFclass.svg` |
| 5I | `Fig5I_EP_by_geneclass_MA.py` | `Fig5I_EP_by_geneclass_MA.svg` |

## Supplementary

| script | output |
|---|---|
| `FigS_saddle.py` | `FigS_saddle_crush.svg`, `FigS_saddle_gc.svg` |
| `FigS_EP_CTCFenhancer_by_GC_MA.py` | `FigS_EP_CTCFenhancer_by_GC_MA.svg` |

## NOTE on 5H — pick one normalisation

5G and 5I are MA-normalised; 5H was specified without. `Fig5H_EP_by_CTCFclass.py` writes
**both** versions so the choice is explicit. Recommend using `Fig5H_EP_by_CTCFclass_MA.svg`
so all three E-P panels share a normalisation — the class ordering is identical either way,
only the offset differs (MA shifts all bars up ~0.010).

## Running

Each panel script is standalone and writes into `svg/`. They set `PROJ=/Zulu/jordan/alcoholATAC`
internally, so they can be run from any directory:

    python3 "Fig5B_compartment_by_GC.py"

They read pre-computed tables in `bensalcohol_out/loops/`. To regenerate those from the raw
Hi-C / ChIP data, run `upstream/` in this order:

    build_anchors.py        -> ctcf_anchor_annot.bed        (CTCF sites + ZF5 triplet + motif score)
    build_looptable.py      -> loop_table.tsv               (SIP loops + O/E + anchor ZF5)
    convergent.py           -> conv_anchor_zf5.tsv          (convergent-orientation anchor calls)
    panel_apa.py            -> apa_loops.pkl                (APA matrices, all loop sets)
    panel_b_profiles.py     -> panelb_profiles.pkl          (CTCF/RAD21/RNAPII anchor profiles)
    panel_rest.py           -> panels.pkl                   (example region, P(s), AA-by-GC)
    saddle2.py              -> saddle_newcrush.pkl          (saddle, crush-rs 5 kb axis)
    ep_annot.py             -> ep_bins.tsv, genes_classed.tsv
    ep_unbiased.py          -> ep_unbiased.tsv              (654,789 annotation-defined E-P pairs)
    ep_assign_enhancer.py   -> ep_unbiased_EP.tsv           (which anchor is the enhancer)
    ma_normalise.py         -> ep_unbiased_MA.tsv           (adds L2E_ma / L2W_ma)

`bullseye.py` is a port of SIPMeta's bullseye transform (Rowley lab `bullseye.py`) used by
5D and 5E; a copy sits both here and in `upstream/`.

## Key methods points for the legend

- Hi-C maps are depth-matched at 592,774,911 contacts (ctrl x0.75, etoh x0.70, T2 native).
- Contact ratios use distance-normalised counts, `(o+1)/(e+1)`.
- E-P pairs are defined by annotation only (enhancer bin x promoter bin, 20 kb-2 Mb).
  They are NOT selected on FitHiC significance: selecting on significance inflates whichever
  condition did the selecting and produced a spurious global up-then-down pattern.
- 5D/5E bullseyes share one colour scale across the three conditions within each panel,
  so the near-identical appearance is the result, not autoscaling.
- 5F permutation null: 20,000 draws, matched on motif score x control CTCF occupancy.
- GC% is used as the duplex guanine measure (each G:C pair contributes one G).
