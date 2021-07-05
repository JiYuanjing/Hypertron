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
				echo "./run_hl_data_v15.csh $a $b $c $d 2 2"
				hadd ht_v15/run18_3gev_ht_hist_pico_merge_centa0_cut2_${a}${b}${c}${d}2_v15.root run18_3gev_ht_hist_pico_rotate_*_centa0_cut2_${a}${b}${c}${d}2_v15.root
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

