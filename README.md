# bicycleModel_osqpMPC
Implementation of the MPC described in [`Combined lateral and longitudinal control for autonomous driving based on Model Predictive Control`](https://webthesis.biblio.polito.it/10667/1/tesi.pdf)

## Dependeces
The project depends only on [`osqp`](http://osqp.readthedocs.io/en/latest/index.html) and [Eigen3](http://eigen.tuxfamily.org/index.php?title=Main_Page).

## Usage
1. Clone the repository

   ```
   git clone https://github.com/robotology/osqp-eigen.git
   ```
2. Build it

   ```
   cd bicycleModel_osqpMPC
   mkdir build
   cd build
   cmake --build .
   ./MPC
   ```
