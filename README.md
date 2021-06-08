# Group Project: Solving the Two Dimensional Poisson Equation

## Introduction

Open Multi-Processing (henceforth, OpenMP) is an application programmming interface (API) capable of supporting shared memory multiprocessing programs. OpenMP makes use of multi-threading, where the program's series of consecutive instructions are broken into discrete sub-threads executed simultaneously on multiple processors. More specifically, OpenMP allows us to *parallelize* the task to ultimately increase the computational efficiency of our program. We demonstrate the utility of parallel programming by numerically solving the two-dimensional Poisson's equation, <img src="/tex/a728fdceeb7edb6d08c57ab25335a6f7.svg?invert_in_darkmode&sanitize=true" align=middle width=133.56565365pt height=26.76175259999998pt/> and comparing the increased computational efficiency to the standard serial programming treatment.

## Theory

As stated, we solve Poisson's equation <img src="/tex/a728fdceeb7edb6d08c57ab25335a6f7.svg?invert_in_darkmode&sanitize=true" align=middle width=133.56565365pt height=26.76175259999998pt/>  where <img src="/tex/f50853d41be7d55874e952eb0d80c53e.svg?invert_in_darkmode&sanitize=true" align=middle width=9.794543549999991pt height=22.831056599999986pt/> is the electrostatic potential and <img src="/tex/6dec54c48a0438a5fcde6053bdb9d712.svg?invert_in_darkmode&sanitize=true" align=middle width=8.49888434999999pt height=14.15524440000002pt/> is the source charge density. In order to solve the problem numerically, we must first introduce a sufficient analytic construction. Ultimately, we will transform the partial differential equation into a matrix problem where the main computational work lies in computing Fourier integrals. 

Since the problem is two-dimensional, we may rewrite Poisson's equation as <p align="center"><img src="/tex/31659c8b0dc6dfd77f3ad90c50fb13cf.svg?invert_in_darkmode&sanitize=true" align=middle width=254.44296074999997pt height=40.11819404999999pt/></p> 

We confine our source charge density to be placed in a metal, rectangular box of side lengths <img src="/tex/9dbcef13f3e6981dfe63f653112a933f.svg?invert_in_darkmode&sanitize=true" align=middle width=18.64161584999999pt height=22.465723500000017pt/> and <img src="/tex/ac6244c7cc0673f3a8d51603f75bcffe.svg?invert_in_darkmode&sanitize=true" align=middle width=18.26684969999999pt height=22.465723500000017pt/>. Placing the walls at <img src="/tex/8436d02a042a1eec745015a5801fc1a0.svg?invert_in_darkmode&sanitize=true" align=middle width=39.53182859999999pt height=21.18721440000001pt/> and <img src="/tex/9dbcef13f3e6981dfe63f653112a933f.svg?invert_in_darkmode&sanitize=true" align=middle width=18.64161584999999pt height=22.465723500000017pt/>, and at <img src="/tex/a42b1c71ca6ab3bfc0e416ac9b587993.svg?invert_in_darkmode&sanitize=true" align=middle width=38.78604674999999pt height=21.18721440000001pt/> and <img src="/tex/ac6244c7cc0673f3a8d51603f75bcffe.svg?invert_in_darkmode&sanitize=true" align=middle width=18.26684969999999pt height=22.465723500000017pt/> imposes the necessary boundary conditions <p align="center"><img src="/tex/121b2aba0948920e666965710a533721.svg?invert_in_darkmode&sanitize=true" align=middle width=164.84245965pt height=92.0098938pt/></p>

For a full analytic treatment, we would normally use separation of variables to obtain two, homogeneous ordinary differential equations (ODEs), apply the boundary conditions, and find the discrete basis functions for each ODE. Putting these solutions together for the non-homogenous equation, we may obtain the expansion of the potential in terms of the basis functions which automatically satisfy the imposed boundary conditions: <p align="center"><img src="/tex/446111e69ab8e69beb7d8ffa09e7eff6.svg?invert_in_darkmode&sanitize=true" align=middle width=368.4128217pt height=49.73538075pt/></p>

Note that these basis functions have been chosen so as to be already normalized. In this basis, the problem simplifies to <p align="center"><img src="/tex/132295f99071c92887c636369deed287.svg?invert_in_darkmode&sanitize=true" align=middle width=261.1920333pt height=49.315569599999996pt/></p>

to be solved for <img src="/tex/3343d1e4776a8c1cc78c3aa3e1c3c557.svg?invert_in_darkmode&sanitize=true" align=middle width=30.808809899999993pt height=14.15524440000002pt/> and where the source term is now

<p align="center"><img src="/tex/4ec3243c4b95787028f2d53b801683ad.svg?invert_in_darkmode&sanitize=true" align=middle width=447.1595766pt height=44.749102199999996pt/></p>

Finally, for the source term we utilize the function <p align="center"><img src="/tex/59c9ff9b238c60d6db45edb314c713c8.svg?invert_in_darkmode&sanitize=true" align=middle width=343.05987539999995pt height=42.07871745pt/></p>
## Methods 

### Serial Version

Here we explain the overall functionality of the program through individual detailing of the specific modules.

#### Modules
1. `main.f90`  
* This module contains the overhead for the program. Here, we call the relevant subroutines to read the namelist file, caclulate the coefficients, write the potential, and print the program run time to the user. Moreover, this module includes a series of logical if statements for running either a series or parallel time test and naming the file for which results are to be written to.
2. `potential.f90` 
* Included functions: `rho_xy`, `rho_mn_int`, `phi_potential`, Included subroutines: `calculate_coefficients`
* This module is responsible for building the solution to Poisson's equation for a charge distribution held within a 2D conducting box. Ultimately, this module builds the electrostatic potential as a function of x and y. To accomplish this, we construct the source term, rho_xy, and then perform a 2D integration over x and y to find rho_mn. From rho_mn, we find the fourier coefficients c_mn. We then sum c_mn multiplied by the basis functions over m,n to construct the potential, phi(x,y). The evaluation of phi(x,y) for all x and y specified within the bounds is computed and recorded to a data file in `read_write.f90`.
3. `quadrature.f90`
* Included functions: `booles_rule`, `booles_quadrature`
* This module contains two functions used to evaluate a volume integral. The numerical integration technique is Boole's quadrature which utilizes 5 equally spaced points to approximate an integral across a given interval in the `booles_rule` function. In `booles_quadrature`, the corresponding integrals for each interval are summed to produce the full integral.
4. `read_write.f90`
* Included subroutines: `read_input`, `write_potential`, `write_wall_clock_time`, `print_wall_clock_time`
* This module contains the subroutine `read_input` intended for reading the provided namelist, as well as establishing default values, in order to run the program with intended parameters. In addition, this module also contains the subroutine `write_potential` for evaluating the solution phi(x,y) at all x and y locations and writing the results to a file. Moreover, this module contains two different routines related to the overall computation time of the program. `print_wall_clock_time` prints the wall clock time to screen upon program completion. `write_wall_clock_time` records the program run time to a file along with the number of threads. In `main`, specific instructions are included to detail how to carry out time tests for parallelization or series execution.

### Parallel Version

The parallelization maintains all the of the same functionality as the serial version. The only major difference is the activation of OpenMP by inclusion of the flag '-fopenmp' in `makefile`, as well as the inclusion of line 207 from module `potential.f90`.

To maximize efficiency, we implement parallelization in the most computationally expensive routines of the program --- two routines total. The first of which is the computation of the coefficients in `calculate_coefficients`. This routine acts as an overhead that ultimately assembles the m x n fourier coefficient matrix. Within this routine, the indices m and n (both up to n_max) are looped over to compute rho_mn and c_mn. The main thread assigns discrete slices of this nested loop to individual subthreads, each of which calls the `rho_mn_int` routine to compute their own double integrals (across all x and y) for a specified section (m and n) of the rho_mn and thus c_mn matrix. The second computationally expensive routine is `write_potential` located in the module `read_write.f90`. Here, the main thread is responsible for assembling the full potential and writing it to a file. It accomplishes this by calling the `phi_potential` function within the module `potential.f90`, where this function assembles the previously computed fourier coefficients and defines the potential as a function of x and y. Thus, the main thread assigns sections of the nested loop (corresponding to x,y locations) to subthreads to compute slices of the potential. The main thread then assembles the full potential as a matrix where i,j correspond to x,y and the value at a given index is the value of the potential at that location.

## Results

### Convergence criteria

One may choose `n_max` based on the level of detailed resolution sought after:

* To accurately illustrate the structure of the potential, the authors recommend an `n_max` between 20 - 40. The overall computation time is on the order of 10^-2 seconds, assuming use of parallelization. Below, we include examples of the 2D and 3D projections of the potential for `n_max = 20`.

<p float="middle">
  <img src="/src/images/2Dnmax20.png" width = "48%"/>
  <img src="/src/images/3Dnmax20.png" width = "48%"/> 
</p>

* For increased resolution, the authors recommend an `n_max` of 70, where increasing values of `n_max` beyond this threshold give little additional resolution. The overall computation time is on the order of 10^-1 seconds. Below, we include examples of the 2D and 3D projections of the potential for `n_max = 70`.

<p float="middle">
  <img src="/src/images/2Dnmax70.png" width = "48%"/>
  <img src="/src/images/3Dnmax70.png" width = 48%"/> 
</p>

## Conclusion

* The program serves as an efficient tool to gain insight into Poisson's equation as well as illustrate the utility of parallelization.

* Overall the program works as intended and is able to perform the calculations of potential in a conducting 2D box with a given charge distribution. Our parallelization is able to improve the speed of the calcualtions as shown in the results section. However, it does seem to level off as we include thread numbers > 3. This could be due to multiple reasons: 
  * The amount of available processing bandwidth 
  * The time testing scripts could have contributed to additional run times
  * Efficiency of the parallelization, meaning there could have been better choices in parallelizations
* For scalability it does not appear there is much to gain using more resources in our parallelization. If we were given more time to assess the problem, we may be able to speed up computation time and possibly find a solution that is scalable. We could do so by utlizing nested parallelization and searching for increased efficiencies in our serial code.

### Instructions for Usage

Within the `box_parameters.namelist` file, specify user chosen initial conditions. If the user chooses not to utilize the namelist file, default values are provided within the `read_write.f90` module for simple illustration. The authors recommend an `n_max` of 70 or over for detailed resolution of the potential. Due to the integration algorithm of Boole's rule, the value of `n_samples` must to be 1mod(4). Moreover, the values of `rho_zero`, `width`, and `center` must all be nonzero and positive. 

#### Implementing the Parallelization

Ensure the flag "-fopenmp" is included within the `makefile`. In addition, ensure that line 207 of module `potential.f90` is *not* commented out. Moreover, ensure all logicals in `main.f90` are set to `.false.`. Assuming all specifications are saved, navigate to your directory and type "make" into terminal. Once all files are compiled type "./poisson box_parameters.namelist" if one wishes to proceed with the provided namelist data, or "./poisson" if one wishes to proceed with default values. Press enter. The program will print to screen the run time as well as the number of threads used. Results are written to a file called `box_parameters.dat`. To visualize the data, open `group_project_analysis.ipynb`. Ensure that this file is located in the same directory as `box_parameters.dat`. At the top of the `.ipynb` file, click run for the first five sections of the program. The first plot illustrates the two dimensional projection of the potential onto the metal conducting box. The second plot illustrated a three dimension representation of the potential, where the z-axis corresponds to the value of the potential and the x and y axis correspond to the length and width of the metal box. Notice that the dark grey plane within the 3D plot illustrates the zero location of the conductor. 

#### Time Testing 

* Parallel Time Testing:

To perform a comparison of the program run time v.s. number of threads, change the logical `run_write_time` in `main.f90` to ".true.". Moreover, ensure that the logical `series` is ".false" and the logical `parallel` is ".true.". After saving the changes, open terminal and navigate to your directory. Utilize the time testing script `time_test_parallel.sh` by typing "sh time_test_parallel.sh" and pressing enter. The shell will execute the program automatically and run the program 10 times for each thread number assignment from 1 - 4. The results will be printed to a file `thread_v_time.dat`. After each trial run, the program will print to screen the run time and the number of threads used. Thus, the program will print out 10 cycles of run times for a thread number of 1, 10 cycles of run times for a thread number of 2, etc. 

* Series Time Testing:

To compare this data with the series program, change the logical `series` to `.true` and the logical `parallel` is `.false.`. Moreover, remove the flag `-fopenmp` in the `makefile` as well as comment out line 270 of module `potential.f90`. After saving the changes, open terminal and navigate to your directory. Now, utilize the time testing script `time_test_serial.sh` by typing "sh time_test_serial.sh" and pressing enter. Again, the shell will execute the program automatically and run the program ten times. However, the automatic thread assignment for the series execution is 0. Thus, the program will print to screen the run time 10 times. This data is saved to a file `serial_time.dat`. 

* Comparing the Two:

To visualize the comparison of the serial v.s. parallel time, open `group_project_analysis.ipynb` and run the last three segments of the program. Below, we include an example output. One notices that although computation time decreases with increasing thread number, the efficiency of the program only minimally increases with thread number assignment beyond 2. Below, we include an example of computation time v.s. thread number for `n_max = 70`. 

<p align="center">
  <img width="400" src="/src/images/threadcomp.png">
</p>



