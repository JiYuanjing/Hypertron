#! /bin/bash
#for data
# root -b -q readtree.C\(\"data_rotate.root\",0,\"fout_H3L_rotate.root\"\,0\)
root -b -q readtree.C\(\"ana_tree.root\",0,\"fout_H3L_data.root\",0,1\)
root -b -q readtree.C\(\"ana_tree.root\",0,\"fout_H3L_data_ME.root\",0,0\)
# for i in `ls list/*.list_*`; do 
#   root -b -q readtree.C\(\"${i}\",0,\"${i}_data.root\"\,0\)
# done
# root -b -q readtree.C\(\"$1\",0,\"$2\"\,0\)
#for mc
# root -b -q readtree.C\(\"Lambda_tree_mc_v2.root\",1,\"fout_Lambda.root\"\)
# root -b -q readtree.C\(\"H3L3b_tree_mc_v2.root\",0,\"fout_H3L.root\"\)
# root -b -q readtree.C\(\"data/H3L3b_tree_mc_quasi.root\",0,\"fout_H3L3b_quasi.root\"\)
# root -b -q readtree.C\(\"data/H3L3b_tree_mc_phase.root\",0,\"fout_H3L3b_phase.root\"\)
# root -b -q readtree.C\(\"data/Lambda_tree_mc_June20.root\",0,\"fout_Lambda_wLaD.root\",-20\)
# root -b -q readtree.C\(\"data/Lambda_tree_mc_June20.root\",1,\"fout_Lambda.root\",1\)
# root -b -q readtree.C\(\"Lambda_tree_mc_rotate_June20.root\",0,\"fout_Lambda_wLaD_rotate.root\",-20\)
