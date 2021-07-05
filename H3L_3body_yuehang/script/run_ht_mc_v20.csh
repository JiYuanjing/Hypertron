#!/bin/csh
set a=0
set b=0
set c=0
set d=0
set i=0
while ( $d < 3 )
	while ( $c < 3 )
		while ( $b < 3 )
			while ( $a < 3 )
				echo "$a $b $c $d 2 2"
				root -b -q process_tree_tmva_ht_kf_vmv20.C+\(0,2,4,0,${a},${b},${c},${d},2\)
				mv run18_3gev_ht_hist_mu_1x_8_newreweight_centa0_wgt4_cut2_${a}${b}${c}${d}2_mve_nohepid_vmv20.root ht_v20
				@ a = $a + 1
				end
			@ a = $a - 3
			@ b = $b + 1
			end
		@ b = $b - 3
		@ c = $c + 1
		end	
	@ c = $c - 3
        @ d = $d + 1
end
@ d = $d - 3

