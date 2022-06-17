# !/bin/bash
# Author: Melyssa Minto


# running homer on zic peaks classified as activators and repressors
# 
./helper_run_homer.sh activating ../../peak_gene/activating/ 500
./helper_run_homer.sh repressive ../../peak_gene/repressive/ 500
./helper_run_homer.sh P7_peaks_DOWNGenes ../../peak_gene/early_activating/ 500
./helper_run_homer.sh P60_peaks_UpGenes ../../peak_gene/late_activating/ 500
./helper_run_homer.sh P60_peaks_DOWNGenes ../../peak_gene/late_repressive/ 500
./helper_run_homer.sh P7_peaks_UpGenes ../../peak_gene/early_repressive/ 500
./helper_run_homer.sh early ../../peak_gene/early/ 500
./helper_run_homer.sh late ../../peak_gene/late/ 500

 ./helper_meme.sh early ../../peak_gene/early/
mv ../../results/invivo/MEME/early/ ../../results/invivo/MEME/early_looped/

./helper_meme.sh late ../../peak_gene/late/
mv ../../results/invivo/MEME/late/ ../../results/invivo/MEME/late_looped/

./helper_meme.sh activating ../../peak_gene/activating/
./helper_meme.sh repressive ../../peak_gene/repressive/
./helper_meme.sh P7_peaks_DOWNGenes ../../peak_gene/early_activating/
./helper_meme.sh P60_peaks_UpGenes ../../peak_gene/late_activating/
./helper_meme.sh P60_peaks_DOWNGenes ../../peak_gene/late_repressive/
./helper_meme.sh P7_peaks_UpGenes ../../peak_gene/early_repressive/
