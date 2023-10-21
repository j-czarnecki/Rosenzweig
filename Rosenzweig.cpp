#include <iostream>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <fstream>
#include <vector>
#include <array>
#include <numeric>

typedef int (*f_equation_system)(double, const double*, double*, void*);
typedef int (*f_jacobian)(double, const double*, double*, double*, void*);

//Global data initialization
const int N_EQUATIONS = 2;
const int N_PARAMS = 5;
const int N_EXTREMA = 100;

// std::vector<std::vector<double>> maxima(N_EQUATIONS);
// std::vector<std::vector<double>> minima(N_EQUATIONS);
std::array<std::array<double, N_EXTREMA>, N_EQUATIONS> maxima;
std::array<std::array<double, N_EXTREMA>, N_EQUATIONS> minima;
std::array<int, N_EQUATIONS> maxima_counter {};
std::array<int, N_EQUATIONS> minima_counter {};


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

bool check_convergence(const double& conv, const double* y, double y_prev[][2]) {
    bool conv_flag = true;
        for (int i = 0; i < N_EQUATIONS; i++){
            if (abs(y[i] - y_prev[i][0]) > conv){
                conv_flag = false;
                break;
            }                
        }
    return conv_flag;
}

void check_extremum(const double* y, double y_prev[][2]){
    for (int i = 0; i < N_EQUATIONS; i++){
        if (y_prev[i][0] > y[i] && y_prev[i][0] > y_prev[i][1]){
            maxima.at(i).at(maxima_counter.at(i)%N_EXTREMA) = y_prev[i][0];
            //std::cout << "Maxima at " << i << "\t" << maxima.at(i).at(maxima_counter.at(i)) << std::endl; 
            maxima_counter.at(i) ++;
        }
        else if(y_prev[i][0] < y[i] && y_prev[i][0] < y_prev[i][1]){
            minima.at(i).at(minima_counter.at(i)%N_EXTREMA) = y_prev[i][0];
            //std::cout << "Minima at " << i << "\t" << minima.at(i).at(minima_counter.at(i)) << std::endl; 
            minima_counter.at(i) ++;
        }
    }
    return;
}

void save_extrema(std::ofstream& OutputExtrema){
    int index {};
    for (int i = 0; i < N_EQUATIONS; i++){
        index = (maxima_counter.at(i) > 100) ? 100 : maxima_counter.at(i);
        OutputExtrema << std::accumulate(maxima.at(i).begin(), maxima.at(i).begin() + index, 0.0)/N_EXTREMA; //maxima mean
        OutputExtrema << "\t";
        index = (minima_counter.at(i) > 100) ? 100 : minima_counter.at(i);
        OutputExtrema << std::accumulate(minima.at(i).begin(), minima.at(i).begin() + index, 0.0)/N_EXTREMA; //minima mean
        OutputExtrema << std::endl;
    }
}

void save_node(std::ofstream& OutputNode, double* y){
    for (int i = 0; i < N_EQUATIONS; i++){
        OutputNode << y[i];
        OutputNode << "\t"; //gurantees NULL at reading second column
        OutputNode << std::endl;
    }
}


void solve_system(double dt, double t_max, double t0, double* y0, double eps_abs, double eps_rel, void* params, \
                const double& conv, std::ofstream& OutputTimeDependence, std::ofstream& OutputCycleOrNode, f_equation_system system, f_jacobian jacobian){
    //GSL procedures allocation
    gsl_odeiv2_system my_equation = {system, jacobian, N_EQUATIONS, params};
    gsl_odeiv2_driver* driver = gsl_odeiv2_driver_alloc_y_new (&my_equation, gsl_odeiv2_step_rk8pd, dt, eps_abs, eps_rel);
    
    //initialize simulation parameters
    double y[N_EQUATIONS] {};
    double y_prev[N_EQUATIONS][2] {};
    double t = t0;
    for (int i = 0; i < N_EQUATIONS; i++){
        y[i] = y0[i];
    }
    int N_steps = int(t_max/dt);
    
    int status {}; //checking for errors

    //convergence checking
    bool conv_flag {}; 

    //Main loop
    for (int i = 0; i < N_steps; i++){
        for (int i = 0; i < N_EQUATIONS; i++){
            y_prev[i][1] = y_prev[i][0];
            y_prev[i][0] = y[i];
        }

        status = gsl_odeiv2_driver_apply (driver,  &t, (i+1)*dt, y); //perform evolution
        
        //std::cout << "Evolution is at time " << t << std::endl;
        OutputTimeDependence << t << "\t" << y[0] << "\t" << y[1] << std::endl;

        if (status != GSL_SUCCESS) {
            printf ("error: driver returned %d at time %f\n", status, t);
            break;
        }

        check_extremum(y, y_prev);

        conv_flag = check_convergence(conv, y, y_prev);
        if (conv_flag){
            std::cout << "Convergence reached" << std::endl;
            break;
        }      

    } //end of main loop

    if (conv_flag){
        save_node(OutputCycleOrNode, y);
    }    
    else{   //if simulation did not converge we can assume that it ended in a cycle
        std::cout << "Convergence not reached, saving cycles" << std::endl;
        std::cout << "Stored maxima: " << maxima.at(0).size() << std::endl;
        save_extrema(OutputCycleOrNode);
    }

    gsl_odeiv2_driver_free (driver);
    return;
}


int main(int argc, char *argv[]){

    try{
        if (argc != 4 + N_EQUATIONS + N_PARAMS){
            throw argc;
        }
    }
    catch(int argc){
        std::cout << "Bad number of params: " << argc << " required: " << 4 + N_EQUATIONS + N_PARAMS << std::endl;
        std::cout << "RUN TERMINATED \n";
        return -1;
    }
    // double dt = 1e-2;
    // double t_max = 1e4;
    double t0 = 0;
    double eps_abs = 1e-8;
    double eps_rel = 1e-8;
    double convergence = 1e-7;

    double y0[N_EQUATIONS] = {0.4, 0.7};
    double params[N_PARAMS] {3. , 1., 1., 2.0, 1.};

    double dt = std::stod(argv[1]);
    double t_max = std::stod(argv[2]);
    for (int i = 0; i < N_EQUATIONS; i++){
        y0[i] = std::stod(argv[3 + i]);
    }

    // params = {p , m , beta, b, alpha_11}
    for (int i = 0; i < N_PARAMS; i++){
        params[i] = std::stod(argv[3 + N_EQUATIONS + i]);
    }
    std::string path = argv[3 + N_EQUATIONS + N_PARAMS];


    std::string filename = path + "/TimeDependence.dat";
    std::ofstream OutputTimeDependence(filename);
    std::string filename_final = path + "/FinalState.dat";
    std::ofstream OutputCycleOrNode(filename_final);
    solve_system(dt, t_max, t0, y0, eps_abs, eps_rel, params, convergence, OutputTimeDependence, OutputCycleOrNode, function, jacobian);

    return 0;
}
