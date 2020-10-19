from random import random,choice
from math import exp
Temp=2.3              # Temperature
n=20                  # Sites per edge for n x n system
n3=n*n*n              # Precalculate n*n*n
nlist=list(range(n))  # List of sites 
ntrials=200000        # Number Trials 
nequil=100000         # Equilibration steps

# Initialize sums for averages
E_sum=0.0
E2_sum=0.0

# Create initial matrix of spins all spin up
spins=[[[1 for i in range(n)] for j in range(n)] for k in range(n)]

# Run simulation
for trial in range(1,(ntrials+nequil+1)):
    # Randomly pick a site
    i=choice(nlist)
    j=choice(nlist)
    k=choice(nlist)

    #Calculate the change in energy if we flip this spin
    deltaE=2.*(spins[i][j][k]*\
               (spins[i][j][(k+1)%n]+spins[i][j][(k-1+n)%n]+\
                spins[i][(j+1)%n][k]+spins[i][(j-1+n)%n][k]+\
                spins[(i+1)%n][j][k]+spins[(i+-1+n)%n][j][k]))

    #Flip the spin using Metropolis MC
    if exp(-deltaE/Temp)>random():   
        spins[i][j][k]=-spins[i][j][k]
    else:
        deltaE=0.0

    # Calculate system energy ONCE 
    if trial == nequil:
        energy=0.0  
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    energy-=(spins[i][j][k]*\
                             (spins[i][j][(k+1)%n]+spins[i][j][(k-1+n)%n]+\
                              spins[i][(j+1)%n][k]+spins[i][(j-1+n)%n][k]+\
                              spins[(i+1)%n][j][k]+spins[(i-1+n)%n][j][k]))
        energy/=n3

    # Update energy based on deltaE per spin
    if trial > nequil:
        energy+=2*deltaE/n3
        E_sum+=energy
        E2_sum+=energy*energy

E_ave=E_sum/ntrials
E2_ave=E2_sum/ntrials
Cv=1./(Temp**2)*(E2_ave-E_ave*E_ave)
print('%8.4f %10.6f %10.6f'%(Temp, E_ave, Cv))
