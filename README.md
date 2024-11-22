This repository contains the codes for the paper "Implicit-explicit time discretization schemes for a class of semilinear wave equations with nonautonomous dampings" by Zhe Jiao, Yaxu Li and Lijing Zhao, cf. https://arxiv.org/abs/2410.20784 


Running the code
1. Make sure that you have installed deal.II (release 9.2.0), cf. https://www.dealii.org/9.2.0
2. open a terminal session in this folder
3. call "cmake -DDEAL_II_DIR=/path/to/deal.II"
4. call "make release" to compile the files. 
5. call "make run" to run the program, or you can put all the files into VS code to run the program.  
6. After that, all data of numerical experiments are generated and saved in the "results" folder.
7. use the Python skript result_plot.py to picture the figures.

Acknowledgements

This project is partially supported by the National Natural Science Foundation of China (No.12272297). Moreover, we have incorporated the spatial discretization part of the codes from the paper "An implicit–explicit time discretization scheme for second-order semilinear wave equations with application to dynamic boundary conditions" by Marlis Hochbruck and Jan Leibold, and we are grateful for their contribution to our project.
