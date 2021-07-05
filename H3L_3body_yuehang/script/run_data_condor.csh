#!/bin/csh
echo 'cut 0: '$1''
echo 'cut 1: '$2''
echo 'cut 2: '$3''
echo 'cut 3: '$4''
set i=0
while ( $i < 8 )
	echo "processing cut:$1,$2,$3,$4, batch:$i"
	root -b -q process_tree_tmva_hl_kf_vmv11.C+\(1,2,0,$i,$1,$2,$3,$4\)
	@ i = $i + 1
end
hadd run18_3gev_hl_hist_pico_merge_loose_8_centa0_cut2_${1}${2}${3}${4}_merge_v11.root run18_3gev_hl_hist_pico_rotate_*_loose_8_centa0_cut2_${1}${2}${3}${4}_v11.root

