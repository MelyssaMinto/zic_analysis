# !/bin/bash
# Author: Melyssa Minto


# running homer and MEME on DE peaks that overlap with DE regions around genes 
# PATH=$PATH:/home/GenomeTools/homer/bin
## These are DA peaks 
# ./helper_run_homer.sh P60vP7_DOWN ../../DiffExp_ZicChIP/ 500
# ./helper_run_homer.sh P60vP7_UP ../../DiffExp_ZicChIP/ 500
# ./helper_run_homer.sh P60vP7_NS ../../DiffExp_ZicChIP/ 500

./helper_meme.sh P60vP7_DOWN ../../DiffExp_ZicChIP/
./helper_meme.sh P60vP7_UP ../../DiffExp_ZicChIP/
./helper_meme.sh P60vP7_NS ../../DiffExp_ZicChIP/