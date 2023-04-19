# Gromacs-ABF

## One of the most important step in gromacs ABF calculation is to define the intermolecular restrain
There are various options available for selecting the atoms for the restraint potential, but the fundamental principle is to opt for atoms that exhibit high rigidity and are less susceptible to relative displacement during the simulation. For small molecules, it is preferable to select atoms that are in proximity to the center of mass position, whereas for proteins, selecting the main chain atoms of amino acids is recommended. In particular, for proteins, the three atoms *A*, *B*, and *C* should ideally have an angle between them ranging from 90 to 120 degrees, with *A* being the closest to the ligand. Additionally, for the ligand atoms *a*, *b*, and *c*, an angle of 90 degrees is preferable if possible, but it is not as strict as the *A*-*B*-*C* angle requirement.

If you are unsure how to select atoms, you can consider choosing the main chain atoms of neighboring amino acids that can form pi-pi interactions. If pi-pi interactions are not possible, atoms near the central part of the small molecule that can form hydrogen bonds could also be considered.

Taking the COVID-19 crystal structure as an example, by Maestro, you can select and define the restrain by yourself following the above general rules

<img width="790" alt="image" src="https://user-images.githubusercontent.com/75652473/232963498-57755fba-87e0-4f9c-a151-5bd329c43f1b.png">

```
Ligand: a:4712; b:4685; c: 4707 
Protein: A:2541; B:2526; C:2525 
r(a-A) = 3.00 
angle baA =159.2 
angle aAB =111.8 
dihedral cbaA =-115.5 
dihedral baAB =87.9 
dihedral aABC =3.9
```
Add the following block to the end of the Gromacs topology file

```
[ intermolecular_interactions]
[ bonds ]
; ai     aj    type   bA      kA     bB      kB
 4712    2541  6      3.00   0.0    3.00   4184.0

[ angles ]
; ai     aj    ak     type    thA      fcA        thB      fcB
 4685   4712   2541   1       159.2     0.0        159.2     41.84
 4712   2541   2526   1       111.8     0.0        111.8     41.84

[ dihedrals ]
; ai     aj    ak    al    type     thA      fcA       thB      fcB
 4707  4685  4712  2541    2       -115.5    0.0     -115.5    41.84
 4685  4712  2541  2526    2        87.9     0.0      87.9     41.84
 4712  2541  2526  2525    2         3.9     0.0      3.9      41.84
```
## Alternatively, it is possible to not define any intermolecular restrain, in this case, you can just skip the restrain correction during later analysis stage, but this is not recommeded.

## Background and introduction

The absolute free energy of binding plays a vital role in protein-small molecule dynamics simulations, particularly in the early stages of drug discovery when the similarity of small molecules is insufficient for calculating relative binding free energies such as FEP. However, the calculation process is complex, requiring a skilled operator, numerous steps, and significant computing power, making it challenging to apply and often costly.

To address these challenges,  this script explicitly designed for Google's Colab computing platform. The script integrates two GROMACS tutorials: Justin's protein complex tutorial and AlchemistryWiki's absolute combined free energy tutorial. The first tutorial is intended for use on a local computer with GROMACS installed, while the second is designed for better computing resources. By simplifying the process and reducing the manual steps for operators, this script provides an opportunity for elementary molecular simulators to perform absolute combined free energy calculations using the free computing power on the network for their research and development activities. This approach can greatly reduce costs and expand the scope of calculating protein-small molecule affinity beyond simulation techniques such as molecular docking.

The script is divided into two parts: the first part is Justin's protein-ligand simulation tutorial, while the second part is the modification inspired by the AlchemistryWiki tutorial. The second part has undergone significant changes to make it suitable for Colab, including changes to the for loop inside run.sh and the lambda folder naming style.  Credits to giribio/MDNotebooks for the software compilation method that was borrowed for this script.

For users who prefer to prepare their ABFE calculations on their laptop, the input files can be uploaded to Colab or AI Studio for use. It is crucial to use the same GROMACS version across different platforms.

ABFE calculations have numerous potential applications, such as evaluating the reliability of docked poses. If a researcher is unsure which of three docked poses is the most likely correct pose, each can undergo ABFE calculations, and the one with the best score can be considered the most probable candidate. For a case study of this approach, see https://pubs.acs.org/doi/10.1021/jacs.6b11467.

## Job control bash script example

## On a 7 cpu cluster with non-mpi version of Gromacs, with GPU
```
# Loop over lambda windows
for (( a=0; a<$n_windows; a++ ))
do
  # Loop over simulations within each window
  for (( b=0; b<$n_sims; b++ ))
  do
    # Create directories for each simulation step
    mkdir -p lambda.$a.$b/ENMIN
    mkdir -p lambda.$a.$b/NVT
    mkdir -p lambda.$a.$b/NPT
    mkdir -p lambda.$a.$b/PROD

    # Energy minimization step
    cd lambda.$a.$b/ENMIN
    gmx grompp -f ../../MDP/ENMIN/enmin.$a.$b.mdp -c ../../solv_ions.gro -p ../../topol.top -n ../../index_jz4.ndx -o enmin.tpr
    gmx mdrun -v -stepout 1000 -s enmin.tpr -deffnm enmin -ntmpi 1 -ntomp 7
    cd ../

    # NVT equilibration step
    cd NVT
    gmx grompp -f ../../MDP/NVT/nvt.$a.$b.mdp -c ../ENMIN/enmin.gro -p ../../topol.top -n ../../index.ndx -o nvt.tpr -r ../../solv_ions.gro
    gmx mdrun -stepout 1000 -s nvt.tpr -deffnm nvt -ntmpi 1 -ntomp 7
    cd ../

    # NPT equilibration step
    cd NPT
    gmx grompp -f ../../MDP/NPT/npt.$a.$b.mdp -c ../NVT/nvt.gro -t ../NVT/nvt.cpt -p ../../topol.top -n ../../index.ndx -o npt.tpr -r ../../solv_ions.gro
    gmx mdrun -stepout 1000 -s npt.tpr -deffnm npt -ntmpi 1 -ntomp 7
    cd ../

    # Production run
    cd PROD
    gmx grompp -f ../../MDP/PROD/prod.$a.$b.mdp -c ../NPT/npt.gro -t ../NPT/npt.cpt -p ../../topol.top -n ../../index.ndx -o prod.tpr
    gmx mdrun -stepout 1000 -s prod.tpr -deffnm prod -dhdl dhdl -ntmpi 1 -ntomp 7
    cd ../../
  done
done

```

## On a 40 cpu cluster with mpi version of Gromacs without GPU (This only applicable to gmx_mpi version, with OpenMPI or interlMPI in the environment)

Create a file called 3.sh (you can name as whatever you like as long as it ends with *.sh)

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
        mpirun -np 1 gmx_mpi mdrun -ntomp 40 -v -stepout 1000 -s enmin.tpr -deffnm enmin
        cd ../
        mkdir NVT
        cd NVT
        gmx_mpi grompp -f ../../MDP/NVT/nvt.$a.$b.mdp -c ../ENMIN/enmin.gro -p ../../topol.top -n ../../index.ndx -o nvt.tpr -r ../../solv_ions.gro
        mpirun -np 1 gmx_mpi mdrun -ntomp 40 -stepout 1000 -s nvt.tpr -deffnm nvt
        cd ../
        mkdir NPT
        cd NPT
        gmx_mpi grompp -f ../../MDP/NPT/npt.$a.$b.mdp -c ../NVT/nvt.gro -t ../NVT/nvt.cpt -p ../../topol.top -n ../../index.ndx -o npt.tpr -r ../../solv_ions.gro
        mpirun -np 1 gmx_mpi mdrun -ntomp 40 -stepout 1000 -s npt.tpr -deffnm npt
        cd ../
        mkdir PROD
        cd PROD
        gmx_mpi grompp -f ../../MDP/PROD/prod.$a.$b.mdp -c ../NPT/npt.gro -t ../NPT/npt.cpt -p ../../topol.top -n ../../index.ndx -o prod.tpr
        mpirun -np 1 gmx_mpi mdrun -ntomp 40 -stepout 1000 -s prod.tpr -deffnm prod -dhdl dhdl
        cd ../../
    done
done

```
## Example of PBS script to use the above 3.sh on a HPC cluster

Create 4.pbs that refers 3.sh to run the simulation

```
#!/usr/bin/csh
#PBS -l select=1:ncpus=40
#PBS -l mem=80gb
#PBS -l walltime=12:00:00
#PBS -o job.o
#PBS -e job.e

cd $PBS_O_WORKDIR

module load intel-mpi/2021.7.1 
module load gromacs/2022.3

sh 3.sh

```
## Submit your job by

```
qsub 4.pbs
```
## Speed and efficiency
For the Justin's tutorial example, 3HTB protein-ligand system, the speed on the 40 CPU only cluster, with Gromacs-mpi/2023, is 25-28 ns/day, for the compelx branch, 50-55ns/day for the ligand branch. One calculation can be done in about 1.5 day without any GPU involved.

## 提高Gromacs的运行效率可以参考中文资料 http://bbs.keinsci.com/thread-13861-1-1.html
