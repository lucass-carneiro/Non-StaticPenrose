#!/bin/bash

spacetime=sks

#files=(sks_breakup_d_6 sks_breakup_d_7 sks_breakup_d_8 sks_breakup_d_9)

#files=(sks_breakup_M2_0 sks_breakup_M2_1 sks_breakup_M2_2 sks_breakup_M2_3 sks_breakup_M2_4 sks_breakup_M2_5)
files=(sks_breakup_M2_6 sks_breakup_M2_7 sks_breakup_M2_8 sks_breakup_M2_9)

for file in ${files[@]}; do
    mkdir $file
    cd $file

    ln -s ../bin/grlensing
    ln -s ../lib/
    ln -s ../../configs/grlensing_config.yaml
    ln -s ../../configs/$spacetime/$file.yaml

    mpirun -n 3 grlensing penrose-breakup $file.yaml > stdout.txt 2> stderr.txt &

    cd ../
done

wait