---
layout: page
mathjax: true
permalink: /Project/
---

## Course Project 
1. [Introduction](#intro)
2. [Deadlines](#deadlines)
3. [Calculations](#calcs)
4. [Analysis](#analysis)<br>
	4a. [LiCoO<sub>2</sub> and Al-doped Data](../CompData)
5. [Final Report](#report)


For the Final Project, you will be studying ethylene carbonate adosrption on Lithium Cobalt Oxide (LiCoO<sub>2</sub>) with different dopant metals. Each group will be assigned a metal dopant to work with (see list of Assigned Projects). The students will work in groups of two or three on the same dopant to perform calculations individually, however, complementing each other. Each pair of students will present their results in class that will be critiqued by the other groups. Finally, each group  will jointly write a final report on the combined data. The due date for the final written report is <font color="red">5/8 at 5:00 PM (hard deadline)</font>.

Please make use of the [Piazza](https://piazza.com/) page for troubleshooting, discussions and for sharing results.

Turn in your final report by emailing a PDF file to:

```
alevoj@seas.upenn.edu, antcurto@seas.upenn.edu
```
<a name='intro'></a>

## Introduction ##

Goal: Understand how different metal dopants affect the structure and reactivity of LiCoO2 specifically through the adsorption of ethylene carbonate (EC).

Plan: Use DFT to study dopant effects on structure and changes it causes to EC adsorption on the surface of LiCoO2. Study the electronic structure of LiCoO2 and doped LiCoO2 to see if any trends exist.

### Motivation ###

LiMO<sub>2</sub> (M= Co, Mn, Ni) materials are the most common cathode materials in batteries today. LiCoO2 is by far the most popular of those materials. A combination of effects has lead to investigations to attempt to improve the stability of LiCoO2 cathode materials or find suitable alternatives. Some problems include cost, reactivity, and structural problems that occur during cycling. Decomposition of the liquid electrolyte is a common cause of decreased battery performance. Ethylene carbonate is a popular component of the elctrolyte solution in Li ion batteries. In order to improve battery performance critical importance in improving battery operation is understanding the electrode – electrolyte interface. For high capacity battery materials, electrolytes decomposed causing the battery to decrease in performance. 

<center><img src="../Images/LiIonBatteries.png" alt="window" style="width: 800px;"/><br>
Schematic of Li Ion Batteries
</center>


AIP Conference Proceedings 1597, 26 (2014); https://doi.org/10.1063/1.4878478 Published Online: 17 February 2015

Ongoing work in the [Cabana Group](https://cabana.chem.uic.edu/) at UIC has shown that some metal dopants, specifically Al in LiCoO<sub>2</sub> can increase the lifetime of the batteries. Well defined nanocrystals of LiCoO<sub>2</sub> can be synthesized and subsequnetly tested for battery performance. The figure below shows a schematic of the nanocrystals synthesized at UIC. With this information we can model the individual facets of the nanocrystal and study the effect on EC adsorption metal dopants will have in each region. 

<center><img src="../Images/nanocrytals.png" alt="window" style="width: 800px;"/><br>
Nanocrystals from the Cabana Group
</center>

Using this information we plan a way to look at specific surface facets and the interactions of these facets with EC. The metals we will look at are Ni, Mg, an Ti. [Data for non-doped LiCoO<sub>2</sub> and Al-doped LiCoO<sub>2</sub> will provided for reference.]((../CompData)

<a name='calcs'></a>

## Calculations ##

To get started we need to get the files needed for the final project. Move into your CBE544 directory and run these commands:

```bash
cd
cd CBE544
wget XXXX
tar -zxvf FinalProject.tar.gz

```
THESE FILES WILL BE UPLOADED BY THE END OF THE DAY ON 4/10. THE GITHUB PAGE WILL BE UPDATE AS SOON AS THESE FILES ARE AVAILABLE.

In the FinalProject directory you should see some trajectory files and a directory called scripts which will be used for this project. In `FinalProject` directory create a `M-surf` and/or `M-subsurf` folder (M is the metal you are assigned, please check [Assignment](https://cbe544.github.io/Project_Assignments/)). For example, if you are assignemnt with Ni and you are running a surface calculation, please run the following command to create the directories: 

```bash
cd FinalProject
mkidr Ni-surf
```

Organization for the project is very important so that the data is accesbile once the class is over. We will structure is to be something like this for the 104 surface

```bash
~/CBE544/FinalProject/104/M-surf/noads/
~/CBE544/FinalProject/104/M-surf/ads1/
~/CBE544/FinalProject/104/M-surf/ads2/
~/CBE544/FinalProject/104/M-surf/ads3/
~/CBE544/FinalProject/104/M-surf/clean/bader
~/CBE544/FinalProject/104/M-surf/ads1/bader
~/CBE544/FinalProject/104/M-surf/ads2/bader
~/CBE544/FinalProject/104/M-surf/ads3/bader
~/CBE544/FinalProject/104/M-sub/noads/
~/CBE544/FinalProject/104/M-sub/ads1/
~/CBE544/FinalProject/104/M-sub/ads2/
~/CBE544/FinalProject/104/M-sub/ads3/
~/CBE544/FinalProject/104/M-sub/clean/bader
~/CBE544/FinalProject/104/M-sub/ads1/bader
~/CBE544/FinalProject/104/M-sub/ads2/bader
~/CBE544/FinalProject/104/M-sub/ads3/bader
```
This is an outline of all of the DFT calculations you will need to do, how to organize the files, and where to run each calculation. For the 001 surface it should be something like this:

```bash
~/CBE544/FinalProject/001/CoTerm/noads/
~/CBE544/FinalProject/001/CoTerm/ads1/
~/CBE544/FinalProject/001/CoTerm/ads2/
~/CBE544/FinalProject/001/CoTerm/ads3/
~/CBE544/FinalProject/001/CoTerm/clean/bader
~/CBE544/FinalProject/001/CoTerm/ads1/bader
~/CBE544/FinalProject/001/CoTerm/ads2/bader
~/CBE544/FinalProject/001/CoTerm/ads3/bader
~/CBE544/FinalProject/001/Literm/noads/
~/CBE544/FinalProject/001/Literm/ads1/
~/CBE544/FinalProject/001/Literm//ads2/
~/CBE544/FinalProject/001/Literm/ads3/
~/CBE544/FinalProject/001/Literm/noads/bader
~/CBE544/FinalProject/001/Literm/ads1/bader
~/CBE544/FinalProject/001/Literm/ads2/bader
~/CBE544/FinalProject/091/Literm/ads3/bader
```

#### A note on Magnetism ###
These calculations are spin polarized and therefore the atoms contain non-zero magnetic moments. In order to view these magnetic moments, which are important for consistent calculations, you can open the final trajectory in ase-gui -> View -> Show Labels -> Magnetic Moments. The problem here is the the ase-gui for ase/3.13.0 does not do a good job with making the magnetic moments readable in the image BUT ase/3.9.1 does. This is further complicated by the fact that we cannot open trajecotry files made in ase/3.13.0 in ase/3.9.1 so we must convert the file format. We can do this using the switch.py script. Specifically, move into the direcotry where you wish to write the final trajectory to an ase/3.9.1 trajectory file and run these commands:

```bash
module load ase/3.9.1
python ~/CBE544/FinalProject/switch.py
```
This will write a file called `mag.traj` which will be openable in ase/3.9.1 and contain the structure and magnetic moments from your caclulation. When presenting structures please include the magnetic moments in your images. This script only works for final trajectories and not for bader charges files.

### Task 1: ### 

Once you have accurately completed HW5 you can continue on to the final project. We will use the 104 surface and 001 surface trajectories provided to you in the FinalProject directory (this is only because there are specific starting magnetic strcutrures that we want otherwise the structures that you made could be used) to place a metal dopant on the surface and subsurface (separately, so two total calculations). The locations are shown here as in a top view of the 104 surface. The simplest way to change an atom to the desired dopant is to use ase-gui, click on the atom to change, Edit (or ctrl+Y), and type in the element you want. Be sure to save this new trejactory because ase-gui does not automatically save any changes you make. 

<center><img src="../Images/dopantlocations.png" alt="window" style="width: 800px;"/><br>
Nanocrystals from the Cabana Group
</center>


Once you have substitued the metal dopant you can use the same relax.py script we used for HW5 to relax this system. Copy over the script submit the job. Record the final energies which can be used to determine if the preferred dopant location is surface or subsurface.

### Task 2: ### 
Using these models from Task 1 we can now adsorb EC to the three locations (per system) shown below:

<center><img src="../Images/Adsorptionlocations104.png" alt="window" style="width: 800px;"/><br>
Locations for Adsorption on the 104 surface of LiCoO<sub>2</sub>
</center>

Refer to the [Adsorption page](../ASE/Adsorption) for instructions on how to add the EC adsorbate. We will use a different script than the relax.py for adsorption calculation. We will use the script called relax-ads.py which you can find in `FinalProject/scripts`. Please be sure to use this script for these calculations or you may expereince converngence issues. 

### Task 3: ### 

Once you have converged the systems with and without EC we will do a [bader charge analysis](http://theory.cm.utexas.edu/henkelman/code/bader/).  Inside the directory where your calculations were run make a new directory called bader (by doing `mkdir bader`). Copy into this directory  fin.traj and vasp-ase.sub. Rename your fin.traj to init.traj (`mv fin.traj init.traj`). Copy from the `FinalProject/scripts` directory the badercharge.py script. It will look like this:

```python
#!/usr/bin/env python
from ase import Atoms, Atom
from ase.calculators.vasp import Vasp
from ase.io import read,write
import numpy as np

p=read('init.traj')
calc = Vasp(prec='accurate',
            encut=520,
            xc='PBE',
            lreal='Auto',
            kpts=[4,4,1],
            nsw = 0,
            ibrion = -1,
            ispin = 2,
            amix_mag = 0.800000,
            bmix = 0.000100,
            bmix_mag= 0.000100,
            amix = 0.20000,
            sigma = 0.05000,
            ediff = 2.00e-04,
            ediffg = -2.00e-02,
            algo ='fast',
            ismear = -5,
            nelm = 250,
            ncore = 16,
            lasph= True,
            ldautype = 2,
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
            lcharg = True,
	    laechg= True,
	    gamma=True,
)
calc.calculation_required = lambda x, y: True
p.set_calculator(calc)
pe=p.get_potential_energy()
#####
ana =  Vasp(restart=True)
pend = ana.get_atoms()

forces=pend.get_forces().ravel()
max_force=max([abs(x) for x in forces])

pe = pend.get_potential_energy()
#mag = pend.get_magnetic_moments()

#pend.set_initial_magnetic_moments(mag)
#print mag
write('fin.traj',pend)
```
This script does a static calculation (nsw=0) of the final trajectory from your previous relaxation and writes the files needed to do a bader charge anaylsis. Use the vasp-ase.sub script to submit the badercharge.py script (`sbatch vasp-ase.sub` with the final line `python badercharge.py`). Once the job has finished you can attach the bader charge to each atom by typing

```python
module load ase/3.9.1
python ~/CBE544/FinalProject/bader_get_charge_vasp
```
This while write a new trajectory file called bader_charge.traj that has attached the bader charge of each atom as a magnetic moment. To see this use ase-gui -> View -> Show Labels -> Magnetic Moments. Analyze how the bader charges differ from each system. It is important to load ase/3.9.1 so that the magnetic moments are readable.

### Task 4: ###

Repeat both Task 1 and Task 2 for the 001 surface. The only difference is instead of surface and subsurface we will use Li-terminated vs CoO termianted. Use the sites clearly depicted below to see where to adosrb the EC in relationship to the metal dopant. The metal dopant in this system should always go in the layer of Co on the surface or closest to it (it will be the layer right under the lithium for the lithium terminated case).

<center><img src="../Images/Adsorptionlocations001-Coterm.png" alt="window" style="width: 800px;"/><br>
</center>
<center><img src="../Images/Adsorptionlocations001-Literm.png" alt="window" style="width: 800px;"/><br>
</center>

Run a bader charge analysis on this system as well. See if there are any clear trends between the two systems through things such as dopant location, charge, etc. Compare your system to the LiCoO<sub>2</sub> and the Al-doped system shown on this page. Look for trends between these systems, your own system, and even those of your classmates (if possible)

<a name='deadlines'></a>

## Deadlines ##
1. HW5 Due: Wed 10 April (Each student)
2. Short update (few slides) on completed calcualtions: Wed 17 April during class (1 per group)
3. Final Presentation: Wed 1 May during class (1 per group)
4. Final Paper: Wed 8 May by 5 PM (1 per group)

## Analysis ##

Do an detailed analysis of your system trying to identifty trends in adsorption due to dopant location, charge, surface, or anything else. Compare your data to that of plain LiCoO<sub>2</sub> and Al-doped LiCoO<sub>2</sub> which can be found [here](/CompData).

### Requirements ###

At a minimum you should accomplish the following:

1. Complete the [HW5](../ASE/Getting_Started).
2. Setup a LiCoO<sub>2</sub> surface (104) and calculate adsorption energies for EC adsorption as three sites for two different metal dopant locations (surface and sub-surface).
3. Do a Bader Charge Analysis on metal doped system and metal doped system w/ EC and compare to the provided systems without and dopant and with an Al dopant.
4. Repeat this process on the 001 facet. Instead of doing is for surface and subsurface we will do this for Li-terminated and CoO<sub>2</sub> terminated. 
5. Analysis
    1. How does the metal dopant affect adsorption vs plain LiCoO<sub>2</sub>? vs Al-doped? How is the adosrption different for surface vs sub-surface?
    2. Do you notice and trends from the bader charge analysis that may contribute to the change in adsorption?
6. Report (3~5 pages maximum)


<a name='report'></a>
### Final Report ###

The final report should be in the form of a 3-5 pages long mini paper including figures and tables. Provide one report for each group. Please be succinct and organize it in the following way:

* Introduction (brief) - don’t write too much
* Calculation details
* Results and discussion including analysis.
* Conclusion (brief)

You are welcome to share data amongst your peers to discuss broader trends. 



