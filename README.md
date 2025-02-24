This repository contains the code that I was using to run simulations of singular SDEs, ran either by Brownian motion of fractional Brownian motion. This is linked to the following paper : https://arxiv.org/abs/2401.09970 

Fractional Brownian motion is simulated by an approximative means with Fourier transform. The goal is to study the zero noise problem in the Peano setting, namely:

dx = b(x) dt + \eps dw^H

where b(x) = A^+ | x|^{\gamma} (x>0) - |x|^{\gamma} (x<0)

The main file is main.cpp. We can pass the parameters:
-eps volatility of the noise
-H Hurst index
-g power of the interaction
-ts timestep
--ips if passed then we run interacting particle system
-A multiplying factor on the positive side in the zero noise problem
-cutoff cutoff in the Fourier simulation of fractional Brownian motion
-fname name of the file in which we save the data

If we pass --ips the script will simulate the interacting particle system, similar to the one in https://arxiv.org/abs/2403.05454  

The output file can be studied with a notebook