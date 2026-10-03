#!/bin/sh
# Reprodukce všech výsledků z PAPER.md. Pořadí je důležité: search.py zapisuje search.npz,
# které další skripty čtou.
set -e
for s in octonion.py labels.py baseline.py symbolic.py convergence.py policy.py \
         search.py percentile.py strict.py ablation.py adaptive.py adaptive2.py free.py clouds.py why_slow.py; do
    echo "============================================================"
    echo "== $s"
    echo "============================================================"
    python3 "$s"
done
