---
layout: page
mathjax: true
permalink: /ASE/Getting_Started/
---

# ASE Tutorials
1. [Introduction to ASE](../)
2. [Getting Started with DFT Calculations](../Getting_Started/)
3. [Adsorption](../Adsorption/)

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

To start this tutorial and the exercises that follow, log on to Stampede2 and download the following:
```bash
wget https://cbe544.github.io/CBE544-2018.github.io/ASE/HW5.tar.gz
tar -zxvf HW5.tar.gz
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

Then, the VASP calculator is set up. All parameters related to the electronic structure calculation are included here. The following example shows typical parameters that we use in the group for MXene calculations.

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

<a name='mxenes'></a>

<a name='lattice-constant-determination'></a>

#### Lattice Constant Determination ####

Find the [`lattice-constant-a.py`](Lattice_Constant.py) script in the `lattice/a` folder. This script calculates the different energies of the system as a function of the lattice constant. Before you run this job, make sure you read the comments within to understand what it does.

```python


Remember to change the script name to lattice-constant-a.py in the `vasp-ase.sub` file! Submit the script by running:

```bash
sbatch vasp-ase.sub
```

To proceed with writing this script, you will be modifying the example script provided here: [ASE-Equation of State](https://wiki.fysik.dtu.dk/ase/tutorials/eos/eos.html). Note that the sample script reads 5 configurations from the trajectory, but we have more configurations than that in our calculations. This script can be run on the login node directly. To execute the script you have written, use the command:


**HW 5:** Plot the energies as listed above, and report the DFT lattice constant.

The two-dimensional bulk modulus B describes the compressibility of a two-dimensional sheet (how difficult it is to stretch or compress). Take the lattice script given before, change the given value to the DFT lattice constant, and change the strain value to run from 0.98 to 1.02, with five steps of 0.1. Fit the energies with a quadratic function. Then, calculuate B via:

$$B=S_{0}\frac{d^{2}E}{dS^{2}}$$

where S is the surface area of the sheet (a variable) and S<sub>0</sub> is the true surface area. Note that in this case we are fitting the surface area S, _not_ the lattice constant! The surface area of the MXenes is given by:

$$S=\frac{\sqrt{3}}{2}a^{2}$$

**HW 5:** Report the two-dimensional bulk modulus of Ti<sub>2</sub>C.

<a name='convergence-with-k-points'></a>

#### Convergence with k-Points ####
Next, we will determine how well-converged the total energy is with respect to the number of k-points in each direction. Modify the [`Lattice_Constant.py`](Lattice_Constant.py) script using the lattice parameter obtained from the previous section. Instead of looping over strain values as above, modify the script to keep the same lattice constant and loop over k-points instead. Try using k = (2,2,1), (4,4,1), (6,6,1), (8,8,1), and (10,10,1), and plot the energy as a function of k-points. Pick one and try to justify why it would be a reasonable choice. The relevant k-points will usually be known, since we have consistent settings that we use throughout the group. In principle, one should always check for convergence when working with a new system. (Side note: Think about why the last k-point is always 1).

**HW 5:** Show the k-point convergence plot, your pick for the k-points, and your rationale.

For the Final Project, we will be using (8x8x1) k-points for the (1x1) MXene surfaces, and (4x4x1) k-points for the (2x2) MXene surfaces.

We have also provided the `Lattice_Resize.py` script that reads in a .traj file, and changes the lattice to the correct version resized with the lattice constant provided. In the script, you need to change the lattice constant you want manually. Unlike the rest of the scripts given to you, this script can be run directly from the command line, without using the submission system, by using:

```bash
python Lattice_Resize.py
```

**Next**: move on to [Adsorption](../Adsorption/) to learn about how to add adsorbates on your surface.
