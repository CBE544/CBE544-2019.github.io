---
layout: page
mathjax: true
permalink: /Project/
---

## Course Project 
1. [Introduction](#intro)
2. [Deadlines](#deadlines)
3. [Calculations](#calcs)
4. [Analysis](#analysis)
5. [Final Report](#report)


For the Final Project, you will be studying ethylene carbonate adosrption on Lithium Cobalt Oxide (LiCoO<sub>2</sub>) with different dopant metals. Each group will be assigned a metal dopant to work with (see list of Assigned Projects). The students will work in groups of two or three on the same MXene to perform calculations individually, however, complementing each other. Each pair of students will present their results in class that will be critiqued by another group of two. Finally, each group  will jointly write a final report on the combined data. The due date for the final written report is <font color="red">5/8 at 5:00 PM (hard deadline)</font>.

Please make use of the [Piazza](https://piazza.com/) page for troubleshooting, discussions and for sharing results.

Turn in your final report by emailing a PDF file to:

```
alevoj@seas.upenn.edu, antcurto@seas.upenn.edu
```

<a name='intro'></a>

<--- ## Introduction ##

The thermo-chemical synthesis of ammonia is accomplished through the [Haber-Bosch process](http://en.wikipedia.org/wiki/Haber_process), where nitrogen gas reacts with hydrogen gas via:

$$
\mathrm{N_2+3H_2\rightarrow 2NH_3}
$$

This process is crucial for the industrial production of fertilizers and chemical feedstocks. Typically, an iron catalyst is used to stabilize the bond-breaking of the N<sub>2</sub> species. The reaction can be separated into elementary reaction steps ([Honkala et. al. (2005)](http://dx.doi.org/10.1126/science.1106435) for more details):

$$
\begin{align}
\mathrm{N_{2\,(g)}} &\rightarrow \mathrm{2N*}\\
\mathrm{H_{2\,(g)}} &\rightarrow \mathrm{2H*}\\
\mathrm{N* + H*} &\rightarrow \mathrm{NH*}\\
\mathrm{NH* + H*} &\rightarrow \mathrm{NH_2*}\\
\mathrm{NH_2* + H*} &\rightarrow \mathrm{NH_3*}\\
\mathrm{NH_3*} &\rightarrow \mathrm{NH_{3\,(g)}}
\end{align}
$$

A free energy diagram is illustrated below:

<center><img src="../Images/N2_path.jpg" alt="N2 path" style="width: 450px;"/>
<br>Ammonia synthesis pathway on a Ru catalyst (<a href="http://dx.doi.org/10.1126/science.1106435">Honkala et. al. (2005)</a>)</center>

Due to the high operating pressures and temperatures required for this reaction, alternative catalysts are still needed for this process. [Medford et. al. (2015)](http://dx.doi.org/10.1016/j.jcat.2014.12.033) have suggested that the linear scaling between the dissociation energy of N<sub>2</sub> and its transition state energy prevents most catalysts from achieving a high rate. Assuming that the bond-breaking of N<sub>2</sub> is rate limiting, then traditional metal catalysts have a transition state that is too high in energy. This is illustrated in the filled contour plot below, where the turnover frequency is plotted as a function of the transition state energy of the first N<sub>2</sub> bond breaking (*E*<sub>N-N</sub>) and the dissociation energy (∆*E*<sub>diss</sub>). A catalyst would need to behave differently from these extended surfaces in order to land in a more active region of the map. 

<center><img src="../Images/N2_volcano.png" alt="N2 volcano" style="width: 400px;"/>
<br>Filled contour plot for the turnover frequencies (Singh et. al. (2016))</center>

We will be exploring 2D MXenes where such configurations might be found. 

Your goals for the project will be to: 
(1) explore this reaction to find unique adsorption configurations and possibly more favorable thermodynamics,
(2) explore if the theormodynamics can be changed with functionalization of the MXenes,
(3) explore relations between the reaction intermediates in the pathway,
(4) compare the chemsitry of the carbide and nitride MXenes.

<a name='deadlines'></a>

## Deadlines ##
1. HW5 Due: Wed 10 April (Each student)
2. Short update (few slides) on completed calcualtions: Wed 17 April during class (per group)
3. Final Presentation: Wed 1 May during class (per group)
4. Final Paper: Wed 8 May by 5 PM (per group)

<a name='calcs'></a>

## Calculations ##

For the Final Project, create a `FinalProj_M2X` folder (M2X is the material you are assigned, please check [Assignment](https://cbe544.github.io/Project_Assignments/)) in your `CBE544` directory. For example, if you are assignemnt with Mo2C, please run the following command to create the directories: 

```bash
cd
cd CBE544
wget CBE544FinalProject
tar -zxvf CBE544FinalProject.tar.gz
cd CBE544FinalProject
```

## Analysis ##

## Final Report ##

The final report should be in the form of a 3-5 pages long mini paper including figures and tables. One report for each group. Please be succinct and organize it in the following way:

* Introduction (brief) - don't write too much
* Calculation details
* Results and discussion
* Conclusion (brief)

You are welcome to share data amongst your peers to discuss broader trends. 

**If you need the energy of the fixed clusters, they are available [here](../Fixed_Lattice_Clusters/energies.txt).**

<a name='grading'></a>

## Grading ##

* 30% exercises
* 20% write-up
* 20% kinetics
* 30% calculations

<a name='reqs'></a>

## Requirements ##

At a minimum you should accomplish the following:

1. Complete the [three exercises](../ASE/).
2. Setup a LiCoO<sub>2</sub> cluster and a (111) surface and calculate adsorption energies for all intermediates.
3. Calculate transition states for the first step N<sub>2</sub> dissociation) using the fixed bond-length method. Extra credit for calculating the hydrogenation barriers.
4. Vibrational frequency and free energy calculations (initial, transition, and final states, and all adsorbed intermediates). 
5. Analysis
    1. Optimal adsorption sites (relation to transition states)
    2. Kinetic rate analysis
6. Report (3~5 pages maximum)

-->
-->
