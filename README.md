# WAND-PIC
![logo](https://github.com/tianhongg/WAND-PIC/blob/master/Resource/Logo.png)


**W**akefield **A**cceleration a**n**d **D**irect-laser acceleration -- **P**article **i**n **C**ell Simulation.

**WAND-PIC** is a 3D parallel quasi-static particle-in-cell simulation code. Developed by Tianhong Wang in 2019, at Dr. Gennady Shvets' Group, Cornell University.

You are welcome to clone or download this repository. :smile: Please also send an email to the author Tianhong Wang(tw474@cornell.edu). We'd like to keep track of user numbers and affiliations.

## Features

* 2D partitioning of the transverse domain.
* Adjustive longitudinal step according to the speed of plasma trajectories.
* Advanced quasi-static equations, there is no need of "predictor-corrector" scheme.
* In-house parallel Multigrid Solver. One for all equations.
* Sub-cycling of macro beam particles, captures sophisticated wave-particle resonant interaction.
* Parallel output and reload.
* Modularity and Scalability.

###  Structure
![logo](https://github.com/tianhongg/WAND-PIC/blob/master/Resource/Code_Stru.png)

## Prerequisites

* A MPI Library
```
impi or mpich
```

* A Parallel NetCDF Library
```
pnetcdf
```

* A Brain. :thought_balloon:

* A Good Weather. :sun_with_face:

* A Good Computer :computer: and Electricity. :zap:


## Compiling
A simple & working Makefile is included [Makefile](Makefile). The 'make' should work on most machines, and your are welcome to write your own one.


###  Tested Environment
**WAND-PIC** has been tested on four different server/supercomputers, from few cores to over 4000 cores, including:

* **TACC STAMPEDE2:** SKX Nodes [Intel Xeon Platinum 8160 ("Skylake")] 
```
intel/18.0.2; impi/18.0.2; pnetcdf/1.11.0
```
* **TACC STAMPEDE2:** KNL Nodes [Intel Xeon Phi 7250 ("Knights Landing")] 
```
intel/18.0.2; impi/18.0.2; pnetcdf/1.11.0
```
* **TACC LONESTAR5:** Nodes [Intel Xeon E5-2690 v3 ("Haswell")] 
```
intel/18.0.2; cray_mpich/7.7.3; pnetcdf/1.8.0
```
* **Small Server:** Nodes [Intel Xeon(R) E7-4830 v4] 
```
intel/18.0.2; impi/18.0.2; pnetcdf/1.11.0
```

(** when running with pnetcdf older than 1.10.0??, parallel output may not work properly. This happens on Lonestar5 occasionally for unknown reasons.)


## Running
Example: running on a server with 36 cores
```
mpiexec -n 36 ./WAND
```


## Developing

I started to work on this project in early 2019, partially for fun, did everything from scratch. And I am currently the only author of this project:scream_cat:. Due to the lack of workforce and limited energy & time, I may improve and update the code slowly in the future:snail::snail::snail::snail:. Hidden bugs are waiting to be found, and many functions need to be added.
This version of WAND-PIC is the first version I finished recently (around July-2019) and it’s already been used in my research. During the development of WAND-PIC, I made many compromises in the code, many places need to be improved. For example, the Multigrid solver I developed, it needs further testing and improving.
Suggestions and feedback are welcome.


## Task Lists
- [x] Periodic Boundary Condition (finished, need test and upload)..
- [ ] More and Better MG Solver Smoother..
- [ ] High-order Pondermotive Potential..
- [ ] OpenMP..
- [ ] Load Balance..


## Author

* **Tianhong Wang (Cornell University)**(tw474@cornell.edu) 



## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* The author of **WAND-PIC** would like to thank all members in Dr. Gennady Shvets' Plasma Group at Cornell University. 
* Thank Dr. Vladimir N. Khudik for deriving related equations.
* Thank Dr. Roopendra Singh Rajawat for inventing the name "WAND-PIC".

#### Reference
[About some of the algorithms](https://aip.scitation.org/doi/abs/10.1063/1.4999629)

[About DLA_1](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.114.184801); [DLA_2](https://aip.scitation.org/doi/abs/10.1063/1.5036967)

[About our group](https://shvets.aep.cornell.edu)

```
||----------------------------------------------------------------------------------||
||----------------------------------------------------------------------------------||
||                                                                                  ||
||               __        ___    _   _ ____        ____ ___ ____                   ||
||               \ \      / / \  | \ | |  _ \      |  _ \_ _/ ___|                  ||
||                \ \ /\ / / _ \ |  \| | | | |_____| |_) | | |                      ||
||                 \ V  V / ___ \| |\  | |_| |_____|  __/| | |___                   ||
||                  \_/\_/_/   \_\_| \_|____/      |_|  |___\____|                  ||
||                                                                                  ||
||----------------------------------------------------------------------------------||
||--  (W)akefield (A)cceleration a(n)d (D)LA - (P)article (i)n (C)ell Simulation  --||
||----------------------------------------------------------------------------------||
||---Author-----------           : Tianhong Wang                --------------------||
||---Starting---------           : Jan-11-2019                  --------------------||
||---Email------------           : tw474@cornell.edu            --------------------||
||---Group------------           : Dr. Gennady Shvets' Group    --------------------||
||---Copyright--------           : (C) 2019 by Tianhong Wang    --------------------||
||----------------------------------------------------------------------------------||
||----------------------------------------------------------------------------------||
```