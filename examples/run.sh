#!/bin/sh

N=3 # (1D dimension, number of spins is N**dim)
dim=2 
J=0.5
h=1
gamma=0.9
T=0.2
nsamples=1000000

python generate_pairs.py $N $dim

python tfim.py $dim $J $h $gamma

make

./main $nsamples $N $dim $J $h $gamma $T
