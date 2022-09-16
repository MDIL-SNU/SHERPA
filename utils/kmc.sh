#!/bin/sh

echo "0" > time.dat
sed -i "s/RESTART .*/RESTART     = 0/g" INPUT
if [ -e "XDATCAR" ]; then
    rm XDATCAR
fi

for i in {1..10000}; do
    # command for SPS
    mpirun -np $nprocs ./SPS

	sed '1,3d' gen_$i/Event.log > log
	reaction=$(awk '{printf( "%s, ", $1 );} END { printf( "\n" ); }' log)
	rate=$(awk '{printf( "%s, ", $3 );} END { printf( "\n" ); }' log)

cat << EOF > ran.py
import random

List = [$reaction]
print(random.choices(List, weights=[$rate])[0])
EOF

    index=$(python ran.py)

cat << EOF > time.py
import sys
import math
import random

print(float(sys.argv[1]) - math.log(random.random())/sum([$rate]))
EOF

    old_time=`tail -n1 time.dat`
    new_time=$(python time.py $old_time)
    echo $new_time >> time.dat

    sed -i "s/INIT_CONFIG.*/INIT_CONFIG \= .\/gen_$i\/Final_$index.POSCAR/g" INPUT
    sed -i "s/OUTPUT_DIR.*/OUTPUT_DIR  = .\/gen_$(($i+1))/g" INPUT
    sed -i "s/RESTART .*/RESTART     = 1/g" INPUT
    sed -i "s/RESTART_DIR.*/RESTART_DIR = .\/gen_$i/g" INPUT

    cat gen_$i\/POSCAR >> XDATCAR

    rm log ran.py time.py
done
