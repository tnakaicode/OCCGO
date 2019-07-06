# OCCGO

Geometrical Optics with Python-OCC

## Branch

* OCC-18.1
  * pythonocc-core==18.1
    * ONLY windows
    * I CANNOT build pythonocc-core==18.2 in Windows
* OCC-18.2
  * pythonocc-core==18.2
    * Build in linux
  * dealii
    * dealii-6
      * step-6
        * Solve the PDE on the current mesh
        * Estimate the error on each cell using some criterion
        * Mark those cells that have large errors for refinement
        * Mark those cells that have particularly small errors for coarsening
        * Refine and coarsen the cells so marked to obtain a new mesh
        * Repeat the steps above until the overall error is sufficiently small
      * setp-6.1 (step-15)
    * dealii-28
    * dealii-54
