# **Implementing Spectral Methods for Partial Differential Equations in Julia**
----------------------------

written by Ernesto Barraza-Valdez with the help of Dr. David Kopriva.

I would like to give a huge thank you to Dr. Kopriva for giving the time, advice, and mentorship to complete these tutorials so that others can more easily access these versatile computational methods. Without Dr. Kopriva's generous guidance, I wouldn't have made it past Chapter 1. 

In the spirit of his compassionate help, please feel free to email or submit issues/discussions/comments.

- Kopriva, D. A. (2009). Implementing spectral methods for partial differential equations: Algorithms for scientists and engineers. Springer Science & Business Media.

-----------------------------

This repository uses Julia to follow Dr. David Kopriva's book: *Implementing Spectral Methods for Partial Differential Equations*

Jupyter notebooks (via Jupyater Lab) are provided for each along with a library that one can use as one goes deeper into the lessons. 

It is split up into Part I, Part IIa, PartIIb, PartIIc, and an extra Part III. All parts are done in 1D for learning purposes. 

* **Part I:** This goes through Chapters 1-3 in the textbook. Sets up the Julia environment with all the packages necessary. The algorithms in these three chapters are developed in Julia

* **Part 2a:** Goes through Chapters 4 in the textbook, Survey of Spectral Approximations. The Fourier, Chebyshev, and Legendre Collocation, Modal Galerkin, Nodal Continuous Galerkin, and Nodal Discontinuous Galerkin Method are reviewed.

* **Part 2b:** Goes through Chapters 5 in the textbook, Spectral Approximations on the Square. The Fourier, Chebyshev, and Legendre Collocation, Modal Galerkin, Nodal Continuous Galerkin, and Nodal Discontinuous Galerkin Method are reviewed in 2D. The wave equation is reviewed.

* **Part 2c:** Goes through Chapters 8 in the textbook, Spectral Elements Method. The Spectral Element Nodal Continuous Galerkin and Nodal Discontinuous Galerkin Methods are reviewed in 1D

* **Part 3:** More 1D examples on the Nodal Discontinuous Galerkin Spectral Element Method. The Euler-Acoustic Equations and Maxwell's Electromagnetic equations for transverse waves are examined. Additionally, a look on how to implement source functions throught the **External State Function** or through Source Functions in the time derivative (Ampere's law current density source)

---------------------

![image](https://github.com/Implementing-Spectral-Methods/Implementing-Spectral-Methods/assets/34816295/dfa53669-7c54-43c9-bdbb-d090f84017f9)


----------------------

# Usage

--------------------
The Jupyter notebooks are standalone and can be used with any Julia version 10+. 

* For Part 1 all package requirements are in the notebook.

* For Part 2a, one should download the notebook **Part 2 Chapter 4** along with the **Part1** submodule

* For Part 2b, one should download the notebook **Part 2 Chapter 5** along with the **Part2a** submodule

* For Part 2c, one should download the notebook **Part 2 Chapter 8** along with the **Part2b** submodule

## Installing Julia

If you're new to Julia, please follow the instructions for installation in the following link:

https://julialang.org/downloads/

## Jupyter Lab with `IJulia.jl`

### Approach 1:

*  To get Jupyter Lab (and access to Jupyter Notebooks via Jupyter Lab), `IJulia.jl` needs to be added to Julia. Follow the instructions of the following link:  https://julialang.github.io/IJulia.jl/stable/manual/installation/

*  Start Julia
*  Type `]` to get into package mode. Then:  `add IJulia`\
  After `IJulia` is installed, use the `backspace` key to get out of package mode. Then type:\
  `using IJulia`\
  If it's the first time, it will compile and build `IJulia`\
  After this use the command:\
  `installkernel("Julia", "--depwarn=no")`\
  This will install Julia for use in Jupyter Lab and notebooks.\
  After the `IJulia` kernel is installed, you can open a notebook or Jupyter Lab. I recommend using Jupyter Lab.
*  To open Jupyter Lab, type into Julia the following: `jupyterlab()`\
*  Jupyter Lab should open in a browser or a link should be outputted that you can paste in your browser.\
  If you are using SSH you will need to add a port and host. You may need to change the `IJulia ` configurations and settings

### Approach 2: Recommended

Julia is great for fast computation but it's still lacking in sophisticated plotting. Having Python to visualize data is still much faster and easier. Thus, the recommendation is to install Anaconda to obtain Python and Jupyter Lab. 

*  Install Anaconda (Miniconda is recommended) follow the instructions here:\
  https://docs.anaconda.com/miniconda/miniconda-install/

*  Once conda is installed you can install Jupyter lab following the links below:\
  https://jupyterlab.readthedocs.io/en/stable/getting_started/installation.html \
  Which basically shows: \
  `conda install -c conda-forge jupyterlab`

*  Once conda and jupyter lab are installed you can start Julia. In the Julia, (following the `IJulia.jl` link) type:\
  `]` to get into package mode. Then:\
  `add IJulia`\
  After `IJulia` is installed, use the `backspace` key to get out of package mode. Then type:\
  `using IJulia`\
  If it's the first time, it will compile and build `IJulia`\
  After this use the command:\
  `installkernel("Julia", "--depwarn=no")`\
  This will install Julia for use in Jupyter Lab and notebooks.\
  After the `IJulia` kernel is installed, you can exit Julia with `exit()`.\
*  With conda installed, you can now open Jupyter Lab using: `Jupyter Lab` command. Depending on your system configurations, a browser window will open with Jupyter Lab or you will need to add some configurations for SSH. 
