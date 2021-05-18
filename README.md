# Kohn, Leibovici & Szkup (2021)

[![Build Status](https://travis-ci.com/vjimenezg/KohnLeiboviciSzkup3.jl.svg?branch=master)](https://travis-ci.com/vjimenezg/KohnLeiboviciSzkup3.jl)
[![Coverage](https://codecov.io/gh/vjimenezg/KohnLeiboviciSzkup3.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/vjimenezg/KohnLeiboviciSzkup3.jl)

This repository contains the complete code used to compute equilibrium steady-state and transition dynamics of the model in [No Credit, No Gain: Trade Liberalization Dynamics, Production Inputs, and Financial Development](https://drive.google.com/file/d/1oDM3Ru-gkF8I4HTdgxnBAiA2BNXvwDfG/view?usp=sharing), by [David Kohn](https://sites.google.com/site/davidkohn16/home), [Fernando Leibovici](https://www.fernandoleibovici.com/) & [Michal Szkup](https://sites.google.com/view/michal-szkup).


The [Online Appendix](https://docs.google.com/viewer?a=v&pid=sites&srcid=ZGVmYXVsdGRvbWFpbnxkYXZpZGtvaG4xNnxneDo3OTBlOTc0ZTExM2VhZTc5) contains all the derivations of the model.


---
### Code directory

* **[/src](/src)** contains all the code used to compute the numerical results presented in the paper, given parameter values (calibration is currently only in MATLAB and available upon request). It comprises the following:

1. **[/src/KLS3.jl](/src/KLS3.jl)** is the main script

2. **[/src/parameters_settings.jl](/src/parameters_settings.jl)** contains the parametrization and settings (including a tauchen function that discretizes the productivity shocks process

3. **[/src/functions.jl](/src/functions.jl)** contains all the functions needed to compute the steady states (KLS3_measure, KLS3_staticproblem, KLS3_dynamicproblem, KLS3_simulate, KL3_GE_par!, KLS3_GE_par)

4. **[/src/functions_trans.jl](/src/functions_trans.jl)** contains all the additional functions needed to compute the transition between two steady states (KLS3_measure_trans, KLS3_staticproblem_period2_altTiming, KLS3_dynamicproblem_trans_vec_t, KLS3_simulate_trans, KLS3_transition_vec2!, KLS3_transition_vec2)

5. **[/src/graphs.jl](/src/graphs.jl)** contains the code to generate the 2D figures in the paper

6. **[/src/welfare.jl](/src/welfare.jl)** contains all the code related to the welfare analysis (**INCOMPLETE**, missing the 3D graphs).


---
### Installation and Use

1. Follow the instructions to [install Julia](https://docs.junolab.org/latest/man/installation/#.-Install-Julia)

2. Download and Install [Atom](https://atom.io/) (you can also use another text editor that supports Julia, but Atom is recommended)

3. In Atom, follow the instructions to [install Juno](https://docs.junolab.org/latest/man/installation/#.-Install-Juno) and similarly, install the language-julia package.

4. Open the [Julia REPL through Atom] (https://docs.junolab.org/latest/man/basic_usage/) and install the package. To do this, enter package mode with `]` (note that this cannot be copy-pasted on OSX; you need to type it to enter the REPL mode) and run

    ```julia
    add https://github.com/vjimenezg/KohnLeiboviciSzkup3.jl.git
    ```

   4a. **Optional, but strongly encouraged**: To install the exact set of packages used here (as opposed to using existing compatible versions on your machine), first run the following in the REPL

      ```julia
      using KohnLeiboviciSzkup3 # will be slow the first time, due to precompilation
      cd(pkgdir(KohnLeiboviciSzkup3))
      ```
      Then, enter the package manager in the REPL (using `]` as before) and type
      ```julia
      ] activate .; instantiate
      ```
      
