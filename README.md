# Gromacs-ABF
Absolute free energy of binding is a very important part of protein-small molecule dynamics simulation. Compared with relative free energy of binding such as FEP, it is more suitable for the relatively more upstream stage of drug discovery, that is, the similarity of small molecules is not enough. Large enough to be able to calculate the relative binding free energy. However, it is similar to FEP and other relative free energy calculations. The calculation process is complicated, requires higher requirements for the operator, many steps, and requires higher computing power, which severely restricts its scope of application, and paid software The price is very expensive, and it cannot be widely used for the time being. I wrote this script specifically for Google’s free computing platform colab. It integrates two existing tutorials based on Gromacs, namely Justin’s protein complex tutorial and AlchemistryWiki’s absolute combined free energy tutorial. The original two tutorials are for the local computer that has Gromacs installed, and the second Absolute Free Energy tutorial is that the user has better computing resources by default. And this script strives to reduce the manual steps of the operator, so that elementary molecular simulators have the opportunity to truly use the free computing power on the network to carry out the method of absolute combined free energy, and apply it in their own scientific research or commercial research and development activities. It can greatly reduce the cost of use, so that the calculation of protein-small molecule affinity is not completely limited to simulation techniques such as molecular docking.

I modified the run.sh in the original script with the overall calculation progress empty, and used a simple and clear for loop statement based on two variables a and b, so that the user can modify ab by himself to control the free energy calculation time not to exceed The 12-hour limit for colab. But it needs to be emphasized that this script is constantly being improved. Therefore, at present, there are still many steps that require the operator to actively participate and intervene. You can use the protein 3htb to go through the whole process first, and then apply it to your own protein-small molecule system. quantaosun@gmail.com

Absolute binding free engergy calculation by Gromacs 2021 GPU. Compilation and Calculations on Google Colab.

Created by Quantao 知乎账号 奋进的涛 https://www.zhihu.com/people/qutesun

This template can be divided into two parts, the fisrt part is the same as Justin's tutorial of protein ligand simulation http://www.mdtutorials.com/gmx/complex/, the second part is my modification inspired by http://www.alchemistry.org/wiki/Absolute_Binding_Free_Energy_-_Gromacs_2016 Note,for the second part, many running commands and overall FEP folder structures have been changed significantly to make it suitable for colab, like the for loop inside run.sh and name style of lambda folders. etc.

Credit to the work of giribio/MDNotebooks as well, since the software compilation method is borrowed there.

If you find installing or running gromacs on online platform is not your preference, you could just use half of this notebook to prepare the ABFE input, then copy the input files to your prefered platform, the only thing to keep in mind is that you need to make sure to use the SAME version of Gromacs across different platforms!

Someone else might prefer to do it the other way round, so you could just prepare your ABFE in your laptop as per the procedures in this notebook and then upload the input files to run the job on colab or AI studio for people in China.

This was added on July 26 2021. There are numberous potential application for ABFE calculation, one of them is to evaluating the reliablity of the docked poses, ie., if you obtained 3 docked poses and not sure which one is more likely the right one, you could subject the 3 to ABFE one by one, the one with the best score could be argued to be the more likely true pose. See https://pubs.acs.org/doi/10.1021/jacs.6b11467 for a case study.

## On a 8 cpu cluster with non-mpi version of Gromacs, with GPU
```
for (( a = 0; a <=2; a++ ))
   do
     for (( b = 0; b <10; b++ ))
          do
          cd lambda.$a.$b
          mkdir ENMIN
          #cd ENMIN
          #gmx grompp -f ../../MDP/ENMIN/enmin.$a.$b.mdp -c ../../solv_ions.gro -p ../../topol.top -n ../../index_jz4.ndx -o enmin.tpr
           #gmx mdrun -v -stepout 1000 -s enmin.tpr -deffnm enmin -ntmpi 1
           #cd ../
           #mkdir NVT
           cd NVT
            gmx grompp -f ../../MDP/NVT/nvt.$a.$b.mdp -c ../ENMIN/enmin.gro -p ../../topol.top -n ../../index.ndx -o nvt.tpr -r ../../solv_ions.gro
           gmx mdrun -stepout 1000 -s nvt.tpr -deffnm nvt
           cd ../
           mkdir NPT
            cd NPT
            gmx grompp -f ../../MDP/NPT/npt.$a.$b.mdp -c ../NVT/nvt.gro -t ../NVT/nvt.cpt -p ../../topol.top -n ../../index.ndx -o npt.tpr -r ../../solv_ions.gro
            gmx mdrun -stepout 1000 -s npt.tpr -deffnm npt -ntmpi 1
             cd ../
            mkdir PROD
            cd PROD 
            gmx grompp -f ../../MDP/PROD/prod.$a.$b.mdp -c ../NPT/npt.gro -t ../NPT/npt.cpt -p ../../topol.top -n ../../index.ndx -o prod.tpr
            gmx mdrun -stepout 1000 -s prod.tpr -deffnm prod -dhdl dhdl -ntmpi 1

            cd ../../
            done
      done
```

## Only a 8 cpu cluster with mpi version of Gromacs without GPU

```
#!/bin/bash

for (( a = 0; a <=2; a++ ))
do
    for (( b = 0; b <10; b++ ))
    do
        cd lambda.$a.$b
        mkdir ENMIN
        cd ENMIN
        gmx_mpi grompp -f ../../MDP/ENMIN/enmin.$a.$b.mdp -c ../../solv_ions.gro -p ../../topol.top -n ../../index_jz4.ndx -o enmin.tpr
        mpirun -np 8 gmx_mpi mdrun -v -stepout 1000 -s enmin.tpr -deffnm enmin
        cd ../
        mkdir NVT
        cd NVT
        gmx_mpi grompp -f ../../MDP/NVT/nvt.$a.$b.mdp -c ../ENMIN/enmin.gro -p ../../topol.top -n ../../index.ndx -o nvt.tpr -r ../../solv_ions.gro
        mpirun -np 8 gmx_mpi mdrun -stepout 1000 -s nvt.tpr -deffnm nvt
        cd ../
        mkdir NPT
        cd NPT
        gmx_mpi grompp -f ../../MDP/NPT/npt.$a.$b.mdp -c ../NVT/nvt.gro -t ../NVT/nvt.cpt -p ../../topol.top -n ../../index.ndx -o npt.tpr -r ../../solv_ions.gro
        mpirun -np 8 gmx_mpi mdrun -stepout 1000 -s npt.tpr -deffnm npt
        cd ../
        mkdir PROD
        cd PROD
        gmx_mpi grompp -f ../../MDP/PROD/prod.$a.$b.mdp -c ../NPT/npt.gro -t ../NPT/npt.cpt -p ../../topol.top -n ../../index.ndx -o prod.tpr
        mpirun -np 8 gmx_mpi mdrun -stepout 1000 -s prod.tpr -deffnm prod -dhdl dhdl
        cd ../../
    done
done

```
