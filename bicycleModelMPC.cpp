#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/MatrixFunctions>
#include <OsqpEigen/OsqpEigen.h>
#include <chrono>

using namespace Eigen;
using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;

struct params  
{  
    double tau = 0.5; // s - drive train time constant
    double C_alpha_f = 19000; // Np/rad - cornering stiffnes front
    double C_alpha_r = 33000; // Np/rad - cornering stiffnes rear
    double m = 1575; // kg
    double L_f = 1.2; // m - CoM to front
    double L_r = 1.6; // m - CoM to rear
    double Iz = 2875; // Nms^2 - yaw moment

    // system size
    int nx = 6;
    int nu = 2;
    int ny = 3;
    int nz = 1;

    // MPC STUFF
    double Ts = 0.1; //s - sampling Time
    double Tl = 1; // s - look-ahead time
    int Np = Tl/Ts;
    int Nc = Np;
    int variables = Nc*nu;

    // Weighting
    double R1 = 1; // weighting: delta
    double R2 = 2; // weighting: ax
    double Q1 = 500; // V-ref weight
    double Q2 = 50; // e1 weight
    double Q3 = 1e3; // e2 weight

    // Constraints
    double axMin = -3; // m/s^2
    double axMax = +3; // m/s^2
    double deltaMin = -50.0/180.0*M_PI; // rad
    double deltaMax = +50.0/180.0*M_PI; // rad
};

struct System {
    MatrixXd A;
    MatrixXd B;
    MatrixXd E;
    MatrixXd C;
};

struct QPmatrizen {
    SparseMatrix<double> hessian;
    VectorXd gradient;
};

struct constraints {
    VectorXd lowerBound;
    VectorXd upperBound;
    SparseMatrix<double> linearMatrix;
};

System setDynamicsMatrices(params &data) {
    System cont;
    cont.A = MatrixXd::Zero(data.nx,data.nx);
    cont.B = MatrixXd::Zero(data.nx,data.nu);
    cont.E = MatrixXd::Zero(data.nx,data.nz);
    cont.C = MatrixXd::Zero(data.ny,data.nx);

    // BUILD SYSTEM MATRIX
    // ---------------------------------------------------------------------------
    cont.A(0,0) = -1/data.tau;
    cont.A(1,0) = 1;
    cont.A(2,2) = -(2*data.C_alpha_f + 2*data.C_alpha_r)/data.m; // missing 1/Vx
    cont.A(2,3) = -(2*data.C_alpha_f*data.L_f - 2*data.C_alpha_r*data.L_r)/data.m; // missing 1/Vx - Vx
    cont.A(3,2) = -(2*data.C_alpha_f*data.L_f - 2*data.C_alpha_r*data.L_r)/data.Iz; // missing 1/Vx
    cont.A(3,3) = -(2.0*data.C_alpha_f*data.L_f*data.L_f + 2.0*data.C_alpha_r*data.L_r*data.L_r)/data.Iz; // missing 1/Vx
    cont.A(4,2) = 1;
    cont.A(4,5) = 1; // missing 1*Vx
    cont.A(5,3) = 1;

    cont.B(0,0) = 1/data.tau;
    cont.B(2,1) = 2*data.C_alpha_f/data.m;
    cont.B(3,1) = 2*data.L_f*data.C_alpha_f/data.Iz;
    
    cont.E(5) = -1; // missing 1*Vx

    cont.C(0,1) = 1;
    cont.C(1,4) = 1;
    cont.C(2,5) = 1;

    return cont;
}

System setDiscreteSystem(System &cont, struct params &data, double Vx) {
    System dis;
    dis.A = cont.A;
    dis.B = cont.B;
    dis.E = cont.E;
    dis.C = cont.C;

    dis.A(2,2) = cont.A(2,2)/Vx;
    dis.A(2,3) = cont.A(2,3)/Vx - Vx;
    dis.A(3,2) = cont.A(3,2)/Vx;
    dis.A(3,3) = cont.A(3,3)/Vx;
    dis.A(4,5) = cont.A(4,5)*Vx;
    dis.E(5) = cont.E(5)*Vx;

    MatrixXd As = MatrixXd::Zero(data.nx + data.nu, data.nx + data.nu); // super A
    As.block(0,0,data.nx,data.nx) = dis.A;
    As.block(0,data.nx,data.nx,data.nu) = dis.B;
    As.block(data.nx,data.nx,data.nu,data.nu) = MatrixXd::Identity(data.nu,data.nu);
    MatrixXd expmAsTs = (As*data.Ts).exp();

    As = MatrixXd::Zero(data.nx + data.nz, data.nx + data.nz); // super A
    As.block(0,0,data.nx,data.nx) = dis.A;
    As.block(0,data.nx,data.nx,data.nz) = dis.E;
    As.block(data.nx,data.nx,data.nz,data.nz) = MatrixXd::Identity(data.nz,data.nz);

    dis.A = expmAsTs.block(0,0,data.nx,data.nx);
    dis.B = expmAsTs.block(0,data.nx,data.nx,data.nu);

    expmAsTs = (As*data.Ts).exp();
    dis.E = expmAsTs.block(0,data.nx,data.nx,data.nz);
    
    return dis;
}

QPmatrizen setHessianGradient(System &dis, params &data, VectorXd &xk, VectorXd &curvature, VectorXd &v_ref) {
    QPmatrizen out;

    // auto t1 = high_resolution_clock::now();

    // start with creating F
    MatrixXd F = MatrixXd::Zero(data.Np*data.ny,data.nx);
    for (size_t i = 0; i < data.Np; i++)
    {
        F.block(data.ny*i,0,data.ny,data.nx) = dis.C*dis.A.pow(i + 1);
    }

    // Make Phi_u
    MatrixXd Phi_u = MatrixXd::Zero(data.Np*data.ny, data.Nc*data.nu);
    MatrixXd firstCol = MatrixXd::Zero(data.Np*data.ny, data.nu);
    firstCol.block(0,0,data.ny,data.nu) = dis.C*dis.A.pow(0)*dis.B;
    for (size_t i = 1; i < data.Np; i++)
    {
        firstCol.block(data.ny*i,0,data.ny,data.nu) = F.block(data.ny*(i - 1),0,data.ny,data.nx)*dis.B;
    }
    for (size_t i = 0; i < data.Nc; i++)
    {
        Phi_u.block(data.ny*i, data.nu*i, data.Np*data.ny - data.ny*i, data.nu) = firstCol.block(0, 0, data.Np*data.ny - data.ny*i, data.nu);
    }

    // // Phi_z
    MatrixXd Phi_z = MatrixXd::Zero(data.Np*data.ny, data.Np*data.nz);
    firstCol = MatrixXd::Zero(data.Np*data.ny, data.nz);
    firstCol.block(0,0,data.ny,data.nz) = dis.C*dis.A.pow(0)*dis.E;
    for (size_t i = 1; i < data.Np; i++)
    {
        firstCol.block(data.ny*i,0,data.ny,data.nz) = F.block(data.ny*(i - 1),0,data.ny,data.nx)*dis.E;
    }
    for (size_t i = 0; i < data.Np; i++)
    {
        Phi_z.block(data.ny*i, data.nz*i, data.Np*data.ny - data.ny*i, data.nz) = firstCol.block(0, 0, data.Np*data.ny - data.ny*i, data.nz);
    }

    // auto t2 = high_resolution_clock::now();
    // duration<double, std::milli> ms_double = (t2 - t1);
    // std::cout << "\nHessian runtime: " << ms_double.count() << "ms";

    // // Make big weighting matrix
    Vector2d R(data.R1,data.R2);
    Vector3d Q(data.Q1,data.Q2,data.Q3);
    MatrixXd bigR = MatrixXd::Zero(data.nu*data.Nc, data.nu*data.Nc);
    MatrixXd bigQ = MatrixXd::Zero(data.ny*data.Np, data.ny*data.Np);
    bigR.diagonal() << R.replicate(data.Nc,1);
    bigQ.diagonal() << Q.replicate(data.Np,1);

    MatrixXd reference = MatrixXd::Zero(data.ny*data.Np,1) ; // 1:3:60 --> v_ref / 2:3:60 --> lateral ref e1 / 3:3:60 --> yaw ref e2
    for (size_t i = 0; i < data.ny*data.Np; i = i + data.ny)
    {
        reference(i) = v_ref(i/data.ny);
    }

    // Build hessian and QP problem
    MatrixXd H = MatrixXd::Zero(data.Nc*data.nu,data.Nc*data.nu);
    H = 2.0*(2.0*Phi_u.adjoint()*bigQ*Phi_u + bigR);
    // H = 0.5*(H + H.adjoint()); // Your Hessian is not symmetric --> make it symmetric!

    // Make vector f
    VectorXd f(data.Nc*data.nu,1);
    f = (-1*0.5*(reference - F*xk - Phi_z*curvature).adjoint()*bigQ*Phi_u).transpose();

    //populate sparse hessian matrix (time intensiv)
    SparseMatrix<double> hessian(data.Np*data.nu,data.Np*data.nu);
    for (size_t i = 0; i < data.Np*data.nu; i++)
    {
        for (size_t ii = 0; ii < data.Np*data.nu; ii++)
        {
            double value = H(i,ii);
            if (value != 0)
            {
                hessian.insert(i,ii) = value;
            }                
        }        
    }

    // assign output
    out.hessian = hessian;
    out.gradient = f;
    return out;
}

constraints setLowerUpperBounds(params &data) {
    constraints out;

    // evaluate the lower and the upper inequality vectors
    VectorXd lowerInequality = VectorXd::Zero(data.Np*data.nu);
    VectorXd upperInequality = VectorXd::Zero(data.Np*data.nu);
    for (size_t i = 0; i < data.Np*data.nu; i += data.nu)
    {
        // constraints for first input (a_x)
        lowerInequality(i) = data.axMin;
        upperInequality(i) = data.axMax;

        // constraints for second input (delta)
        lowerInequality(i+1) = data.deltaMin;
        upperInequality(i+1) = data.deltaMax;
    }
    
    // populate linear constraint matrix
    SparseMatrix<double> A(data.Np*data.nu,data.Np*data.nu);
    for (size_t i = 0; i < data.Np*data.nu; i++)
    {
        A.insert(i,i) = -1;
    }
    
    // asssert output
    out.linearMatrix = A;
    out.lowerBound = lowerInequality;
    out.upperBound = upperInequality;
    return out;
}


int main()
{   
    // get params
    params data;

    // OP inputs
    // ---------------------------------------------------------------------------
    VectorXd v_ref(data.Np,1);
    v_ref = VectorXd::Ones(data.Np)*2.7778; // m/s

    VectorXd xk(data.nx,1);
    xk << 0, 0, 0, 0, 0, 0.6568;

    VectorXd curvature(data.Np,1);
    curvature << 0.0530,  0.1152,    0.1264,    0.1385,    0.1514,    0.1650,    0.1794,    0.1942,    0.2093,    0.2244;
    // ---------------------------------------------------------------------------
    
    // System-Matrix
    System cont = setDynamicsMatrices(data);

    // Discrtize system !!
    System dis = setDiscreteSystem(cont,data,0.1);

    // Build Hessian, f, constraint matrix etc.
    QPmatrizen qp_matrizen = setHessianGradient(dis,data,xk,curvature,v_ref);

    // Make constraints
    constraints cons = setLowerUpperBounds(data);

    // OSQP - solve the problem
    OsqpEigen::Solver solver;
    solver.settings()->setWarmStart(true);
    solver.settings()->setVerbosity(true); // disable solver feeback

    solver.data()->setNumberOfVariables(data.Nc*data.nu);
    solver.data()->setNumberOfConstraints(data.Np*data.nu);
    if(!solver.data()->setHessianMatrix(qp_matrizen.hessian)) return 1;
    if(!solver.data()->setGradient(qp_matrizen.gradient)) return 1;
    if(!solver.data()->setLinearConstraintsMatrix(cons.linearMatrix)) return 1;
    if(!solver.data()->setLowerBound(cons.lowerBound)) return 1;
    if(!solver.data()->setUpperBound(cons.upperBound)) return 1;

    // instantiate the solver
    if(!solver.initSolver()) return 1;

    // controller input and QPSolution vector
    VectorXd QPSolution;

    // solve the QP problem
    if(!solver.solve()) return 1;

    // get the controller input
    QPSolution = solver.getSolution();

    std::cout << QPSolution << std::endl;
    
    auto t1 = high_resolution_clock::now();
    // Loop ten times and calcute average time
    for (size_t i = 0; i < 10; i++)
    {
        // Update hessian --> next step
        dis = setDiscreteSystem(cont,data,0.1); // Update first the discrete system ---> velocity update!
        qp_matrizen = setHessianGradient(dis,data,xk,curvature,v_ref); // Build hessian etc
        if(!solver.updateHessianMatrix(qp_matrizen.hessian)) return 1;
        if(!solver.updateGradient(qp_matrizen.gradient)) return 1;

        // Solve again
        if(!solver.solve()) return 1;
        QPSolution = solver.getSolution();

    }

    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms_double = (t2 - t1)/10;
    std::cout << "\nAverage runtime: " << ms_double.count() << "ms\n";
}    

