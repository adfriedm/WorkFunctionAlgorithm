# Work Function Algorithm

## Brief Description
This repository includes a somewhat naive approach to the Work Function
Algorithm for finite general psuedo metric spaces in Python. The implementation
is based on the dynamic programming algorithm given in "The k-server
problem"-Koutsoupias which is a summary paper which includes the algorithm that
Koutsoupias and Papadimitriou introduced earlier.

## Included Files
wfa.py - the class implementation

wfa\_example.py - an example of the algorithm running on a Euclidean grid with
the usual metric. Running will produce PNGs in a figs/ folder.

make\_mpg.sh - generates a MOV file from the figures. Requires ImageMagick and
the ffmpeg plugin.

# Background
## K-Server Problem
Given k servers on a pseudometric space, determine a sequence of configurations
in an online manner so that a given sequence of requests is served by a server
in each corresponding configuration.

## Competitiveness
The algorithm has been proven to be (2k-1)-competitive on general metric spaces
and has been conjectured to be k-competitive. In a practical sense, this means
that no matter how poor a choice requests are, the cost of the Work Function
Algorithm is guaranteed to perform no worse than (2k-1) times worse than an
optimal algorithm that knows all requests ahead of time. It is known that no
deterministic online algorithm can do better than k times the optimal offline.

## Work Function Algorithm
The Work Function Algorithm balances the greedy algorithm (assigning the
closest server to a request) and the retrospective algorithm (moves servers the
optimal configuration assuming the current request is the terminating request).
It is the best proven competitive algorithm known for the k-server problem.
There are flow implementations that are polynomial in k and the number of
points in the space, however this is the original exponential time algorithm
with a number of amoritisations included.

# Future Updates
There will be two main updates in the coming months,

1. Fast implementation based on "A fast implementation of the optimal offline
algorithm for solving the k-server problem"- Rudec et al.

2. Limited memory implementation based on "Work function algorithm can forget
history without losing competitiveness"- Colussi

Both of these will make using the algorithm more practical or even just
feasible for large spaces, sequence lengths and numbers of servers.
