# TopologySensitivity
Implementation of rank-1 inversion-free step to get topology sensitivity information for combinatorial topology optimization


<img src="topoSens-minQ-results.png" height="400" class="center"/>

<em>An example of the structure, current and cost function resulting from evaluation of START.m script.</em>

## Implementation notes

The local step from [1] is fully implemented, including additions and removals. The topology sensitivity evaluation is fully vectorized, starts here from a random seed, and converges to local minima in Greedy sense. Either removals or additions can be disabled. The optimization can be restricted to a part of the structure or its boundary. Two fitness functions are shown: minimization of Q factor, w/o preference on self-resonance. They are presented as the first example in [2]. The fitness function is evaluated with radiation matrix factorized whenever useful. The code works with arbitrary excitation.

## Example

The example utilizes pre-calculated data from method of moments simulation to reduce both the code complexity and the computational cost of the demonstration. The code is fully compatible with the outputs of AToM package [3]. Pre-calculated data are provided for perfectly conducting rectangular region of aspect ratio 2:1, with discrete delta gap feeder placed in the middle. The plate is discretized into 512 triangles and covered with 744 basis functions, 743 of which serves as degrees of freedom for combinatorial topology optimization.

## Initiation and start

Optimization parameters can be set at the beginning of START.m script, which serves as a starting script and runs automatically after pressing "F5". No extra code is required, except of visualisation of current density, which requires AToM [3] installed.

It is advantageous to normalize fitness function traces with fundamental bounds. They can be evaluated with in-house code "FunBo" (Fundamental Bounds Package) [4], which is an add-on to AToM, and can be freely downloaded.


## References

[1] Capek, M., Jelinek, L., Kadlec, P., Gustafsson, M.: Optimal Inverse Design Based on Memetic Algorithm -- Part 1: Theory and Implementation, arXiv preprint, arXiv: xxxxx, pp. 1-12, 2022.

[2] Capek, M., Jelinek, L., Kadlec, P., Gustafsson, M.: Optimal Inverse Design Based on Memetic Algorithm -- Part 2: Examples and Properties, arXiv preprint, arXiv: xxxxx, pp. 1-13, 2022.

[3] Antenna Toolbox for MATLAB (AToM), [on-line]: www.antennatoolbox.com, (2022)

[4] Fundamental Bounds Package (FunBo, AToM add-on), [on-line]: http://antennatoolbox.com/atom#addons (2022)
