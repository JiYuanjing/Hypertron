#!/bin/csh
set a=0
set b=0
set c=0
set d=0
set i=0
while ( $d < 1 )
	while ( $c < 1 )
		while ( $b < 1 )
			while ( $a < 1 )
				while ( $i < 8 )
					echo "processing cut:$a,$b,$c,$d, batch:$i"
					root -b -q process_tree_tmva_hl_kf_vmv12.C+\(1,2,0,$i,$a,$b,$c,$d\)
					@ i = $i + 1
					end
					hadd run18_3gev_hl_hist_pico_merge_loose_8_centa0_cut2_${a}${b}${c}${d}_merge_v12.root run18_3gev_hl_hist_pico_rotate_*_loose_8_centa0_cut2_${a}${b}${c}${d}_v11.root
				@ i = $i - 8
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

