#include <iostream>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <fstream>

#define N_EQUATIONS 2
#define N_PARAMS 5


int function(double t, const double y[], double dydt[], void* params){
    (void)(t);
    double* params_arr = (double*) params;
    double p = params_arr[0];
    double m = params_arr[1];
    double beta = params_arr[2];
    double b = params_arr[3];
    double alpha_11 = params_arr[4];

    dydt[0] = y[0]*(p*beta*y[1]/(1 + b*(beta*y[1])) - m); // consumer
    dydt[1] = y[1]*(1 - alpha_11*y[1]) - y[0]*p*y[1]*beta/(1 + b*beta*y[1]);//resource 1
    return GSL_SUCCESS;
}

int jacobian(double t, const double y[], double* dfdy, double dfdt[], void* params){
    (void)(t);
    double* params_arr = (double*) params;
    double p = params_arr[0];
    double m = params_arr[1];
    double beta = params_arr[2];
    double b = params_arr[3];
    double alpha_11 = params_arr[4];

    gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, N_EQUATIONS, N_EQUATIONS);
    gsl_matrix* jacobian_matrix = &dfdy_mat.matrix;
    gsl_matrix_set (jacobian_matrix, 0, 0, p*beta*y[1]/(1 + b*beta*y[1]) - m);
    gsl_matrix_set (jacobian_matrix, 0, 1, p*y[0]*(beta*(1 + b*beta*y[1]) - beta*beta*b*y[1])/pow(1 + b*beta*y[1],2));
    gsl_matrix_set (jacobian_matrix, 1, 0, p*beta*y[1]/(1 + b*beta*y[1]));
    gsl_matrix_set (jacobian_matrix, 1, 1, 1 - 2*alpha_11*y[1] - p*y[0]*(beta*(1 + b*beta*y[1]) - beta*beta*b*y[1])/pow(1 + b*beta*y[1],2));
    dfdt[0] = 0.;
    dfdt[1] = 0.;
    return GSL_SUCCESS;
}


void solve_system(double dt, double t_max, double t0, double* y0, double eps_abs, double eps_rel, void* params, \
                double conv){


    std::string filename = "Plot_phase.dat";
    std::string output_data {};
    std::ofstream OutputFile(filename);

    gsl_odeiv2_system my_equation = {function, jacobian, N_EQUATIONS, params};
    gsl_odeiv2_driver* driver = gsl_odeiv2_driver_alloc_y_new (&my_equation, gsl_odeiv2_step_rk8pd, dt, eps_abs, eps_rel);
    
    double y[N_EQUATIONS] {};
    double y_prev[N_EQUATIONS] {};
    double t = t0;
    for (int i = 0; i < N_EQUATIONS; i++){
        y[i] = y0[i];
    }

    int N_steps = int(t_max/dt);
    int status {};
    bool conv_flag {};

    //Main loop
    for (int i = 0; i < N_steps; i++){
        for (int i = 0; i < N_EQUATIONS; i++){
            y_prev[i] = y[i];
        }

        status = gsl_odeiv2_driver_apply (driver,  &t, (i+1)*dt, y); //perform evolution
        
        //std::cout << "Evolution is at time " << t << std::endl;
        OutputFile << t << "\t" << y[0] << "\t" << y[1] << "\t" << std::endl;

        if (status != GSL_SUCCESS) {
            printf ("error: driver returned %d at time %f\n", status, t);
            break;
        }
        conv_flag = true;
        for (int i = 0; i < N_EQUATIONS; i++){
            if (abs(y[i] - y_prev[i] > conv)){
                conv_flag = false;
                break;
            }                
        }

        // if (conv_flag){
        //     std::cout << "Convergence reached" << std::endl;
        //     break;
        // }      

    }

    return;
}


int main(){

    double dt = 1e-2;
    double t_max = 1e4;
    double t0 = 0;
    double y0[N_EQUATIONS] = {0.4, 0.7};
    double eps_abs = 1e-8;
    double eps_rel = 1e-8;

    // params = {p , m , beta, b, alpha_11}
    double params[N_PARAMS] {3. , .75, 0.5, 2.0, 0.95};

    double convergence = 1e-12;

    solve_system(dt, t_max, t0, y0, eps_abs, eps_rel, params, convergence);

    return 0;
}
