# Network Creation with Homophilic Agents
We study Network Creation Games with multiple types of homophilic agents and non-uniform edge cost, introducing two models focusing on the perception of same-type
and different-type neighboring agents, respectively.

This repository contains a C++ implementation of the game-theoretic model defined in "Network Creation with Homophilic Agents" paper (see  for the full version or proceedings
of IJCAI-ECAI 22 (the 31st International Joint Conference on Artificial Intelligence and the 25th European Conference on Artificial Intelligence)).

### Input
All parameters of the game are defined in main() (see main.cpp). 

### Output 
For each input of the game, the following files will be generated:
  1. GE\_\<number of nodes>\_n\_\<alpha>\_alpha\_\<number of black nodes>\_of\_black\_\<version of the game>\_from\_\<initial network>\_rand\_\<best or improving responce\>\_\<number of run\>\_run\_global\_\<avg clustering coefficient of the resulting network\>\_avg\_clust\_\<diameter\>\_D\_\<number of steps requiered to converge\>\_steps.dot
    Contains the resulting stable network
  2. \<number of nodes\>\_n\_\<alpha\>\_alpha\_\<number of black nodes>\_of\_black\_\<version of the game\>\_from\_\<initial network>\_rand\_\<best or improving responce\>\_avg\_segregation.dot 
    Average clustering coefficient of each pairwise stable network generated by the algorithm at each run
  3. \<number of nodes\>\_n\_\<alpha\>\_alpha\_\<number of black nodes>\_of\_black\_\<version of the game\>\_from\_\<initial network>\_rand\_\<best or improving responce>\_\<number of run\>\_run\_\<network global segregation\>\_segregation.dot
    Local segregation of each node in a pairwise stable network generated by the algorithm
  4. \<number of nodes\>\_n\_\<alpha\>\_alpha\_\<number of black nodes>\_of\_black\_\<version of the game\>\_from\_\<initial network>\_rand\_\<best or improving responce>\_\<number of run\>\_run\_moves\_timeline.txt 
    A timeline of the algorithm steps
  5. \<number of nodes\>\_n\_\<alpha\>\_alpha\_\<number of black nodes>\_of\_black\_\<version of the game\>\_from\_\<initial network>\_rand\_\<best or improving responce>\_\<number of run\>\_run\_avg\_segr\_timeline.txt 
    A timeline of the network segregation at each step of the algorithm

### Further notes
The current version of the implementation doesn't support both version of the game simultiniously.
To switch between the game versions (DEI- and ICF-NCG), please uncomment the correct line in SchellingNCG::cost_of_the_neighborhood function.

This code was written for a more general setting. Hence, you can see the classes' hierarchy and some extra features that you probably don't need for computing a pairwise stable network.
