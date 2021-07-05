#!/bin/csh
echo 'cut 0: '$1''
echo 'cut 1: '$2''
echo 'cut 2: '$3''
echo 'cut 3: '$4''
echo 'cut 4: '$5''
set i=0
while ( $i < 8 )
	echo "processing cut:$1,$2,$3,$4,$5 batch:$i"
	root -b -q process_tree_tmva_ht_kf_vmv12.C+\(1,2,0,$i,$1,$2,$3,$4,$5\)
	@ i = $i + 1
end

