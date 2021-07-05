#!/bin/csh
echo 'cut 0: '$1''
echo 'cut 1: '$2''
echo 'cut 2: '$3''
echo 'cut 3: '$4''
echo 'cut 4: '$5''
echo 'cut 5: '$6''
set i=0
while ( $i < 8 )
	echo "processing cut:$1,$2,$3,$4,$5,$6 batch:$i"
	root -b -q process_tree_tmva_hl_kf_vmv15.C+\(1,$6,0,$i,$1,$2,$3,$4,$5\)
	@ i = $i + 1
end

mv run18_3gev_hl_hist_pico_centa0_cut${6}_${1}${2}${3}${4}${5}_v15.root hl_v15
hadd hl_v15/run18_3gev_hl_hist_pico_merge_centa0_cut${6}_${1}${2}${3}${4}${5}_v15.root run18_3gev_hl_hist_pico_rotate_*_centa0_cut${6}_${1}${2}${3}${4}${5}_v15.root    

