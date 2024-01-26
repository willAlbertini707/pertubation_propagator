# Orbital Propagator with Pertubations  
William J. Albertini  
Aerospace Engineering Department California Polytechnic San Luis Obispo  

This project is a performative orbital propagator built  
for the AER0 452-01-2338 final project. It takes into account  
orbital pertubations to model orbital changes/decay.  

## To Use  
The propagator program should be compiled using the  
rust package manager with the following command:  

> cargo build --release

The Python wrapper class in env_builder.py is then  
used to load in environment variables and run the  
program.  

The output is a .txt with spacecraft state data. Flags  
can be used to turn on and off pertubations.


