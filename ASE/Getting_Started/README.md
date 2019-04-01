---
layout: page
mathjax: true
permalink: /ASE/Getting_Started/
---

# ASE Tutorials
1. [Introduction to ASE](../)
2. [Getting Started with DFT Calculations](../Getting_Started/)

____

## Getting Started with DFT Calculations ##

In the first exercise, we will be studying lithium cobalt and how to determine their lattice constants, followed by surface relaxation of the (104) surface. For Homework 5, everyone will be studying the same system (104) LiCoO<sub>2</sub>. For the Final Project, you will use the same system but with multiple facets (104 and 001) to study ethylene carbonate adosrption.

## Contents ##

1. [A Typical ASE Script](#a-typical-ase-script)
2. [Lattice Constant Determination](#lattice-constant-determination)
3. [Convergence with k-points](#convergence-with-k-points)
4. [Optimization](#optimization)


<a name='a-typical-ase-script'></a>

### A Typical ASE Script ###

ASE scripts can be run directly in the terminal (in the login node) or submitting to external nodes. Generally, you will be submitting jobs to external nodes and only small scripts will be run on the login node. By default, all output from any submitted script will be written *from the directory where the submission command was executed*, so make sure you are inside the calculation folder before running the submission command.

To start this tutorial and the exercises that follow, log on to Chestnut and download the following:
```bash
wget https://cbe544.github.io/CBE544-2019.github.io/HW.tar.gz
tar -zxvf HW.tar.gz
cd HW5
```

There are two files that are necessary to run jobs on the Chestnut cluster. The first is `vasp-ase.sub`; this is the file that tells the scheduler how much time the job is allowed, how many processors it requires, and other pertinent information. First, notice the comments in the beginning. These include information such as how much time to allocate, the number of nodes required, what the names of the output and error files are, what the name of the job should be, and what your email is. 

```bash
#!/bin/bash

#SBATCH -x node63,node64,node81                 #node to exclude due to problems. Do not change
#SBATCH -p p_alevoj                             #partition to run on. do not change
#SBATCH -N  2 #number of nodes
#SBATCH --tasks-per-node=32                     #do not change
#SBATCH -t 48:30:00 #time limit
#SBATCH -J JOBNAME #job name                    #Name your job here
#SBATCH -o out.%j #screen output                #output file name. %j is job number
#SBATCH -e err.%j #errinfo                      #err file name
#SBATCH --mail-user=EMAL@seas.upenn.edu         #add you email address here to be alerted when job ends
#SBATCH --mail-type=end #notify when job finishes #mail when job ends


#export VASP_GAMMA=true
module load ase-vasp/run        #load vasp and ase

python script-name.py   #name of script to run
```

Finally, the last line ```python script-name.py``` picks the script you want to run. Therefore, you need to change the name of the file depending on which script you are running. We will be using this script later in this section for performing calculations to compute the lattice constant of bulk LiCoO<sub>2</sub>.


Let's look at how a typical ASE script is written. Open the [`relax.py`](energy.py) script. We import all the relevant ASE modules in for this calculation

```python
from ase import Atoms, Atom
from ase.calculators.vasp import Vasp
from ase.io import read,write
import numpy as np
```

`from ase.calculators import Vasp` imports the VASP calculator for the ASE interface, and `from ase.io import read, write` imports the read and write commands for trajectory files.

An existing trajectory can be read in:

```python
# read in the slab
slab = read('LiCoO2.traj')
```

Then, the VASP calculator is set up. All parameters related to the electronic structure calculation are included here. The following example shows typical parameters that we use in the group for LiCoO<sub>2</sub> calculations.

```python
calc = Vasp(prec='normal',	#scf accuracy
            encut=520,		#plane-wave cutoff
            xc='PBE',		#functional
            lreal='Auto',	#sampling space
            kpts=[4,4,1],	#kpoint sampling
            nsw = 99,		#max number of ionic steps
            ibrion = 2,		#ion iteration steps
            ispin = 2,		#spin polarized
            amix_mag = 0.800000,#mixing parameters
            bmix = 0.000100,
            bmix_mag= 0.000100,
            amix = 0.20000,
            sigma = 0.05000,	#smearing
            ediff = 2.00e-04,	#energy difference for scf convergence
            ediffg = -2.00e-02,	#force  cutoff for overall convergence
            algo ='fast',
            ismear = -5,	#smearing type
            nelm = 250,		#max number of electronic steps
            ncore = 16,
            lasph= True,
            ldautype = 2,	#Use Hubbard U
            lmaxmix = 4,
            lorbit = 11,
            ldau = True,
            ldauprint = 2,
            ldau_luj={'Co':{'L':2, 'U':3.32, 'J':0},
                      'Li':{'L':-1, 'U':0.0, 'J':0.0},
                      'O':{'L':-1, 'U':0.0, 'J':0.0}
                      },
            lvtot = False,
            lwave = False,
            lcharg = False,
	    gamma=True,		#center at gamma point
)
```

Finally, the VASP calculator is attached to the `slab` Atoms object, the energy calculation is ran, and the total energy of the system is output in the log file (defined in the `spede_esp.sub` file above). 

Once the scripts and atoms object is set up you can submit a job, using:

```bash
sbatch vasp-ase.sub
```

<a name='lattice-constant-determination'></a>

#### Lattice Constant Determination ####

Find the [`lattice-constant-a.py`](Lattice_Constant.py) script in the `lattice/a` folder. This script calculates the different energies of the system as a function of the lattice constant. Before you run this job, make sure you read the comments within to understand what it does.

The following lines have been added the the beginning of the script to vary to lattice size of the bulk LiCoO<sub>2</sub>. LiCoO<sub>2</sub> is symmetric in two directions therefore our a and b lattice constants are going to the same and both mus tbe changed at the same time. However the c lattice vector can be studied independent of a and b. In addition the unit cell is created on an angle. This can be seen in the b and c initial lattice parameters. 

```python
eps=0.03
a0=2.835
c0=4.71
a=[a0,0,0]
b=[a0/2,(a0/2)*np.sqrt(3),0]
c=[a0/2,a0/(2*np.sqrt(3)),c0]
for X in np.linspace(1-eps,1+eps,7):
	p=read('init.traj')
 	p.set_cell([[i*X for i in a],[j*X for j in b],c],scale_atoms=True)
```

Remember to change the script name to lattice-constant-a.py in the `vasp-ase.sub` file! Submit the script by running:

```bash
sbatch vasp-ase.sub
```

To proceed with writing this script, you will be modifying the example script provided here: [ASE-Equation of State](https://wiki.fysik.dtu.dk/ase/tutorials/eos/eos.html). Note that the sample script reads 5 configurations from the trajectory, but we have 7 in our calculations.  

This script can be run on the login node directly. To execute the script you have written, use the command:

```python
python EOS-script.py
```

This will plot the volumes vs energies and print the volume that related to the minimum energy. The output plot (EOS.png) should show the fitted energies as a function of the volume, with the volume corresponding to the minimum and the bulk modulus displayed on the top. To get the a<sub>DFT</sub> lattice constant take this volume and use this equation:

a<sub>DFT</sub> = {2*Volume}/{(4.71)*sqrt(3)}

Repeat this process with the c lattice constant by going to the lattice/c directory and running the script there. The process for the EOS will be the same but to get the c lattice constant you must use this formula.

c<sub>DFT</sub> = sqrt[(Volume/6.959)^2+1.4175^2+0.818^2]

**HW 5:** Show your Python scripts for the EOS, Plot the Equation of State fits, and report the DFT lattice constants.

<a name='convergence-with-k-points'></a>

#### Convergence with k-Points ####
Next, we will determine how well-converged the total energy is with respect to the number of k-points in each direction. You will be running the kptconv.py script in the k-points folder. Look through the script to understand what its doing. Run this script by submitting a job to an external node as discussed previously. Remember to change the name of the script to execute, in the spede_esp.sub file. Upon completion, the script outputs a convergence plot and prints the total energies as a function of the k-points used in the calculation.

**HW 5:** Show the k-point convergence plot, your pick for the k-points, and your rationale.

#### Optimization ####
Finally, you will be performing a geometry optimization on the 104 surface of LiCoO<sub>2</sub>. To proceed with this exercise, first take a look at the starting structure `LiCoO2-104.traj` in the `relax` folder by using the GUI. You should see a 1x4x6 surface of LiCoO<sub>2</sub>, with the bottom three layers fixed to the bulk positions. Next, take a look at the `relax.py` script discussed previously. You will be using this script for running the surface optimization calculations. Submit the calcualtion using vasp-ase.sub and be sure to change the file name accordingly.

To check on the calculation while it is running use either 

```bash
more out.XXXXXX
````
or 
```bash
less out.XXXXX
```
Once the calculation is completed and you have recieved an email go to the relax directory. The final energy will be printed in the last line of the out.XXXXX file. To see this line use the command `tail out.XXXXX`. Tail will print out the last few lines of the file. You should see something like this:

```bash
Energy = -514.7679737
max_forces = 0.015439
```
Please be sure that the force is below the criteria set. The Energy listed is the energy of the system.

**HW 5:** Report the converged energy of the optimized structure. 

