# quantum-walk
python code for a particle hopping on a one-dimensioanl lattice accroding to a quantum rule for moving left/right. The particle interacts with a memory qubit at the site to determine its next move. This qubit has a memory of previous visits to the site.
The parameters at the top of the file when changed will generate various walk patterns. It will also plot the probability distributions from the interference pattern as 
the particle's wave fucntion spreads across the lattice. 

The walk uses scipy package, and uses optimized routines for matrix handling. The size of these arrays doubles each time a lattice point is added. 

WARNING: Exercise care when using n>11 (number of lattice points) as the matrix size will be huge.

Run it! 
