#! /bin/bash
#for data
# root -b -q readtree.C\(\"data_rotate.root\",0,\"fout_H3L_rotate.root\"\,0\)
# input="ana_tree.root"
# input="data_test3.root"
# input="dataME_test2.root"
# input="test.list"
input="quasiMC2.root"
# root -q -b readMc.C\(\"$input\",0,\"fMC_H3L_0080.root\",1\)
# root -b -q readtree.C\(\"$input\",0,\"fout_H3L_MC_0080_010pt.root\",1,0,0,8\)
# root -q -b readMc.C\(\"$input\",0,\"fMC_H3L_0050.root\",1,3,8\)
root -b -q readtree.C\(\"$input\",0,\"fout_H3L_MC_0050_015pt.root\",1,0,3,8\)
# root -b -q readtreesys.C\(\"$input\",0,\"fout_H3L_MC_0050_015pt_sys.root\",1,0,3,8\)
# root -b -q readtree.C\(\"$input\",0,\"fout_H3L_MC_0050_tight.root\",1,0,3,8\)
# root -b -q readtree.C\(\"$input\",0,\"fout_H3L_MC_0050_010pt.root\",1,0,3,8\)

input="phaseMC.root"
# root -b -q readtreesys.C\(\"$input\",0,\"fout_H3L_phaseMC_0050_015pt_sys.root\",1,0,3,8\)
# root -b -q readtree.C\(\"$input\",0,\"fout_H3L_phase_MC.root\",1,0\)
# input="LambdaMC.root"
input="LambdaMC_mix_pid3.root"
# root -b -q readtree.C\(\"$input\",0,\"fout_Lambda_MC_Cuts_0050_015pt.root\",-20,0,3,8\)
# root -b -q readtreesys.C\(\"$input\",0,\"fout_Lambda_MC_Cuts_0050_015pt_sys_pid2.root\",-20,0,3,8\)
# root -b -q readtreesys.C\(\"$input\",0,\"fout_Lambda_MC_Cuts_0050_015pt_sys.root\",-20,0,3,8\)
# root -b -q readtreesys.C\(\"$input\",0,\"fout_Lambda_MC_Cuts_0050_015pt_sys_pid3.root\",-20,0,3,8\)
# root -b -q readtreesys.C\(\"$input\",0,\"fout_Lambda_MC_Cuts_0050_015pt_sys_mixpid3.root\",-20,1,3,8\)
# root -b -q readtree.C\(\"$input\",0,\"fout_Lambda_MC_Cuts_0050_tight.root\",-20,0,3,8\)
# root -b -q readtree.C\(\"$input\",0,\"fout_Lambda_MC_Cuts_0050_010pt.root\",-20,0,3,8\)
# root -b -q readtree.C\(\"$input\",0,\"fout_H3L_data.root\",0,0\)
# root -b -q readtree.C\(\"$input\",0,\"fout_H3L_data_ME.root\",0,1\)
# root -b -q readtree.C\(\"$input\",0,\"fout_H3L_data_SE.root\",0,2\)
# root -b -q readtree.C\(\"$input\",0,\"fout_H3L_data_check.root\",0,0\)
# root -b -q readtree.C\(\"$input\",0,\"fout_H3L_data_ME_check.root\",0,1\)
# root -b -q readtree.C\(\"$input\",0,\"fout_H3L_data_SE_check.root\",0,2\)

# input="testME.root"
# root -b -q readtree.C\(\"$input\",0,\"fout_ME_test.root\",0,1\)
# root -b -q readtree.C\(\"$input\",0,\"fout_KF_test.root\",0,0\)

# root -b -q readtree.C\(\"$input\",0,\"fout_H3L_data_pcut.root\",0,0\)
# root -b -q readtree.C\(\"$input\",0,\"fout_H3L_data_ME_pcut.root\",0,1\)
# root -b -q readtree.C\(\"$input\",0,\"fout_H3L_data_SE_pcut.root\",0,2\)


# root -b -q readtree.C\(\"ana_tree_3012.root\",0,\"test.root\",0,-1\)


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
