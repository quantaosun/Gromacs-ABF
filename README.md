# Gromacs-ABF

## Update: A potential automatic tool to generate Gromacs intermolecular restrain

https://chemrxiv.org/engage/api-gateway/chemrxiv/assets/orp/resource/item/61dea70f81f3fee545a4825a/original/evaluating-the-use-of-absolute-binding-free-energy-in-the-fragment-optimization-process.pdf (page 8/9)

## One of the most important step in gromacs ABF calculation is to define the intermolecular restrain
There are various options available for selecting the atoms for the restraint potential, but the fundamental principle is to opt for atoms that exhibit high rigidity and are less susceptible to relative displacement during the simulation. For small molecules, it is preferable to select atoms that are in proximity to the center of mass position, whereas for proteins, selecting the main chain atoms of amino acids is recommended. In particular, for proteins, the three atoms **A**, **B**, and **C** should ideally have an angle between them ranging from 90 to 120 degrees, with **A** being the closest to the ligand. Additionally, for the ligand atoms **a**, **b**, and **c**, an angle of 90 degrees is preferable if possible, but it is not as strict as the **A**-**B**-**C** angle requirement.

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
 4712    2541  6      0.30   0.0    0.30   4184.0

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

## After the simulaiton is done, and analysis is finised, we correct the result by the following script, it is critical to pay attention to the r0's is nm not angstrom

```
#!/usr/bin/env python

import math
import sys

#===================================================================================================
# INPUTS
#===================================================================================================

K = 8.314472*0.001  # Gas constant in kJ/mol/K
V = 1.66            # standard volume in nm^3

T      = 300.0      # Temperature in Kelvin
r0     = 0.654      # Distance in nm
thA    = 88.8       # Angle in degrees
thB    = 32.9       # Angle in degrees

K_r    = 4184.0     # force constant for distance (kJ/mol/nm^2)
K_thA  = 41.84      # force constant for angle (kJ/mol/rad^2)
K_thB  = 41.84      # force constant for angle (kJ/mol/rad^2)
K_phiA = 41.84      # force constant for dihedral (kJ/mol/rad^2)
K_phiB = 41.84      # force constant for dihedral (kJ/mol/rad^2)
K_phiC = 41.84      # force constant for dihedral (kJ/mol/rad^2)

#===================================================================================================
# BORESCH FORMULA
#===================================================================================================

thA = math.radians(thA)  # convert angle from degrees to radians --> math.sin() wants radians
thB = math.radians(thB)  # convert angle from degrees to radians --> math.sin() wants radians

arg =(
    (8.0 * math.pi**2.0 * V) / (r0**2.0 * math.sin(thA) * math.sin(thB))
    *
    (
        ( (K_r * K_thA * K_thB * K_phiA * K_phiB * K_phiC)**0.5 ) / ( (2.0 * math.pi * K * T)**(3.0) )
    )
)

dG = - K * T * math.log(arg)

print "dG_off = %8.3f kcal/mol" %(dG/4.184)
print "dG_on  = %8.3f kcal/mol" %(-dG/4.184)

```

## Background and introduction

The absolute free energy of binding plays a vital role in protein-small molecule dynamics simulations, particularly in the early stages of drug discovery when the similarity of small molecules is insufficient for calculating relative binding free energies such as FEP. However, the calculation process is complex, requiring a skilled operator, numerous steps, and significant computing power, making it challenging to apply and often costly.

This script integrates two GROMACS tutorials: Justin's protein complex tutorial and AlchemistryWiki's absolute combined free energy tutorial. 

This workflow is divided into two parts: the first part is Justin's protein-ligand simulation tutorial, while the second part is the modification inspired by the AlchemistryWiki tutorial. For the first part, you need do the Justin procedure twice, one for complex, one for ligand only. **Alternatively, you can use charmm gui as updated inside the charmmgui folder of this repository.**

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

## Example of analysis of 3HTB solvation 

```
Temperature: 300 K

Detailed results in kT (see help for explanation):

 lam_A  lam_B      DG   +/-     s_A   +/-     s_B   +/-   stdev   +/- 
     0      1    0.00  0.00    0.00  0.00    0.00  0.00    0.00  0.00
     1      2    0.00  0.00    0.00  0.00    0.00  0.00    0.00  0.00
     2      3    0.00  0.00    0.00  0.00    0.00  0.00    0.00  0.00
     3      4    0.00  0.00    0.00  0.00    0.00  0.00    0.00  0.00
     4      5    0.00  0.00    0.00  0.00    0.00  0.00    0.00  0.00
     5      6    0.00  0.00    0.00  0.00    0.00  0.00    0.00  0.00
     6      7    0.00  0.00    0.00  0.00    0.00  0.00    0.00  0.00
     7      8    0.00  0.00    0.00  0.00    0.00  0.00    0.00  0.00
     8      9    0.00  0.00    0.00  0.00    0.00  0.00    0.00  0.00
     9     10    0.00  0.00    0.00  0.00    0.00  0.00    0.00  0.00
    10     11    3.42  0.03    0.89  0.02    0.81  0.01    1.33  0.01
    11     12    1.99  0.02    0.62  0.01    0.55  0.00    1.09  0.01
    12     13    1.01  0.02    0.42  0.02    0.38  0.02    0.91  0.01
    13     14    0.29  0.02    0.34  0.01    0.32  0.01    0.82  0.01
    14     15    0.92  0.00    0.08  0.00    0.09  0.00    0.42  0.00
    15     16    0.89  0.00    0.10  0.00    0.11  0.00    0.43  0.00
    16     17    1.67  0.01    0.36  0.01    0.53  0.01    0.91  0.01
    17     18    1.52  0.00    0.42  0.01    0.63  0.02    0.99  0.01
    18     19    1.29  0.01    0.52  0.02    0.81  0.03    1.12  0.01
    19     20    0.96  0.02    0.64  0.01    1.06  0.01    1.28  0.01
    20     21    0.34  0.03    0.96  0.03    1.56  0.06    1.80  0.02
    21     22   -0.18  0.03    0.35  0.02    0.52  0.22    1.21  0.06
    22     29   -3.86  0.91    9.64  0.67   29.54  0.68   45.28  2.60


Final results in kJ/mol:

point      0 -      1,   DG  0.00 +/-  0.00
point      1 -      2,   DG  0.00 +/-  0.00
point      2 -      3,   DG  0.00 +/-  0.00
point      3 -      4,   DG  0.00 +/-  0.00
point      4 -      5,   DG  0.00 +/-  0.00
point      5 -      6,   DG  0.00 +/-  0.00
point      6 -      7,   DG  0.00 +/-  0.00
point      7 -      8,   DG  0.00 +/-  0.00
point      8 -      9,   DG  0.00 +/-  0.00
point      9 -     10,   DG  0.00 +/-  0.00
point     10 -     11,   DG  8.54 +/-  0.07
point     11 -     12,   DG  4.96 +/-  0.05
point     12 -     13,   DG  2.53 +/-  0.04
point     13 -     14,   DG  0.73 +/-  0.04
point     14 -     15,   DG  2.29 +/-  0.01
point     15 -     16,   DG  2.22 +/-  0.01
point     16 -     17,   DG  4.16 +/-  0.03
point     17 -     18,   DG  3.78 +/-  0.01
point     18 -     19,   DG  3.22 +/-  0.03
point     19 -     20,   DG  2.39 +/-  0.06
point     20 -     21,   DG  0.84 +/-  0.08
point     21 -     22,   DG -0.45 +/-  0.07
point     22 -     29,   DG -9.62 +/-  2.28

total      0 -     29,   DG 25.59 +/-  2.25
```
dG = -25.59 KJ/mol (-6.1 Kcal/mol)

