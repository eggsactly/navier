# About

This project provides a cool way to visualize fluid flow around objects. 
The sample program is a slightly modified version of a Navier-Stokes fluid 
dynamics code by Griebel et al. [1].

It's also extended off the repo at 
[github.com/dorchard/navier](https://github.com/dorchard/navier).

# History 
2012/13 - Accelerate (Haskell) and Ypnos (Haskell) implementations in progress.
2020 - Implementation of a bitmap upload interface for airfoil shapes. 

# Building and Running

```
make
./navier
^C
```
^C when you're done running.

Then make the video by:

```
./build-vid output
```

# References 
[1]: Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer, Numerical Simulation in Fluid Dynamics, SIAM, 1998. http://people.sc.fsu.edu/~jburkardt/cpp_src/nast2d/nast2d.html
