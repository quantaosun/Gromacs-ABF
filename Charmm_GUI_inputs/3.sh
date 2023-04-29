#!/bin/bash

for (( a = 0; a <=2; a++ ))
do
    for (( b = 0; b <10; b++ ))
    do
        cd lambda.$a.$b
        mkdir ENMIN
        cd ENMIN
        gmx_mpi grompp -f ../../MDP/ENMIN/enmin.$a.$b.mdp -c ../../step3_input.gro -p ../../topol.top -n ../../index.ndx -o enmin.tpr -r ../../step3_input.gro -maxwarn 7
        mpirun -np 1 gmx_mpi mdrun -ntomp 40 -v -stepout 1000 -s enmin.tpr -deffnm enmin
        cd ../
        mkdir NVT
        cd NVT
        gmx_mpi grompp -f ../../MDP/NVT/nvt.$a.$b.mdp -c ../ENMIN/enmin.gro -p ../../topol.top -n ../../index.ndx -o nvt.tpr -r ../ENMIN/enmin.gro -maxwarn 7
        mpirun -np 1 gmx_mpi mdrun -ntomp 40 -stepout 1000 -s nvt.tpr -deffnm nvt
        cd ../
        mkdir NPT
        cd NPT
        gmx_mpi grompp -f ../../MDP/NPT/npt.$a.$b.mdp -c ../NVT/nvt.gro -t ../NVT/nvt.cpt -p ../../topol.top -n ../../index.ndx -o npt.tpr -r ../NVT/nvt.gro -maxwarn 7
        mpirun -np 1 gmx_mpi mdrun -ntomp 40 -stepout 1000 -s npt.tpr -deffnm npt
        cd ../
        mkdir PROD
        cd PROD
        gmx_mpi grompp -f ../../MDP/PROD/prod.$a.$b.mdp -c ../NPT/npt.gro -t ../NPT/npt.cpt -p ../../topol.top -n ../../index.ndx -o prod.tpr -r ../NPT/npt.gro  -maxwarn 7
        mpirun -np 1 gmx_mpi mdrun -ntomp 40 -stepout 1000 -s prod.tpr -deffnm prod -dhdl dhdl
        cd ../../
    done
done

