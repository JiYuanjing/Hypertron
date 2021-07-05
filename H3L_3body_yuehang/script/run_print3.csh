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
				hadd run18_3gev_ht_hist_pico_merge_loose_8_centa0_cut2_${a}${b}${c}${d}_merge_v11.root run18_3gev_ht_hist_pico_rotate_*_bdt_10_centa0_cut2_${a}${b}${c}${d}_v11.root 
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

