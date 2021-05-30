# bicycleModel_osqpMPC
Implementation of the MPC described in [`Combined lateral and longitudinal control for autonomous driving based on Model Predictive Control`](https://webthesis.biblio.polito.it/10667/1/tesi.pdf) with the [`osqp`](http://osqp.readthedocs.io/en/latest/index.html) solver, [`Eigen3`](http://eigen.tuxfamily.org/index.php?title=Main_Page) and the wrapper [`osqp-eigen`](https://github.com/robotology/osqp-eigen).

## Dependeces
- [`osqp`](http://osqp.readthedocs.io/en/latest/index.html)
- [`Eigen3`](http://eigen.tuxfamily.org/index.php?title=Main_Page)
- [`osqp-eigen`](https://github.com/robotology/osqp-eigen)

## Usage
1. Clone the repository

   ```
   git clone https://github.com/phmandl/bicycleModel_osqpMPC
   ```
2. Build it

   ```
   cd bicycleModel_osqpMPC
   mkdir build
   cd build
   cmake --build .
   ./bicycleModel_osqpMPC
   ```
