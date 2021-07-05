#!/bin/csh
set a=0
set b=0
set c=0
set d=0
while ( $d < 3 )
	while ( $c < 3 )
		while ( $b < 3 )
			while ( $a < 3 )
				echo "processing cut:$a,$b,$c,$d"
				root -b -q process_tree_tmva_hl_kf_vmv13.C+\(0,2,4,0,$a,$b,$c,$d,4\)
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

