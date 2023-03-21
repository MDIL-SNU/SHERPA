# setting #
KMC_STEP=1000
TEMP=300
ATT_FREQ=1e13
INIT_CONFIG=./POSCAR

rm KMC.XDATCAR time.log
sed -i "s|INIT_CONFIG.*|INIT_CONFIG\t= ${INIT_CONFIG}|" INPUT
for ((i = 0; i < ${KMC_STEP}; ++i)); do
    sed -i "s|OUTPUT_DIR.*|OUTPUT_DIR\t= ${i}|" INPUT
    cat $(grep 'INIT_CONFIG' INPUT | awk {'print $3'}) >> KMC.XDATCAR
    mpirun ./SPS_LMP
    ./KMC ${i} $((${i}+1)) ${TEMP} ${ATT_FREQ}
    sed -i "s|INIT_CONFIG.*|INIT_CONFIG\t= ./tmp_POSCAR|" INPUT
    sed -i "s|RANDOM_SEED.*|RANDOM_SEED\t= -1|" INPUT
done
