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
const int N_CONSUMERS = 3;
const int N_RESOURCES = 3;
const int N_EQUATIONS = N_CONSUMERS + N_RESOURCES;
const int N_PARAMS = 3 + N_RESOURCES*N_RESOURCES + N_CONSUMERS*N_RESOURCES;
const int N_EXTREMA = 100;
const double BETA_CONVERGENCE_TIME= 5e4;

double beta_convergence_counter = 0;
// std::vector<std::vector<double>> maxima(N_EQUATIONS);
// std::vector<std::vector<double>> minima(N_EQUATIONS);
std::array<std::array<double, N_EXTREMA>, N_EQUATIONS> maxima;
std::array<std::array<double, N_EXTREMA>, N_EQUATIONS> minima;
std::array<std::array<double, N_EXTREMA>, N_EQUATIONS> maxima_time;
std::array<std::array<double, N_EXTREMA>, N_EQUATIONS> minima_time;
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

int function_double_resource_adaptive(double t, const double y[], double dydt[], void* params){

    (void)(t);
    double* params_arr = (double*) params;
    double p = params_arr[0];
    double m = params_arr[1];
    double b = params_arr[2];
    double alpha_11 = params_arr[3];
    double alpha_12 = params_arr[4];
    double alpha_21 = params_arr[5];
    double alpha_22 = params_arr[6];
    double v = params_arr[7];

    dydt[0] = p*y[0]*(y[3]*y[1] + (1-y[3])*y[2])/(1 + b*(y[3]*y[1] + (1 - y[3])*y[2])) - m*y[0]; //consumer
    dydt[1] = y[1]*(1 - alpha_11*y[1] - alpha_12*y[2]) - (p*y[1]*y[3]*y[0])/(1 + b*(y[3]*y[1] + (1 - y[3])*y[2])); //resource 1 
    dydt[2] = y[2]*(1 - alpha_21*y[1] - alpha_22*y[2]) - (p*y[2]*(1 - y[3])*y[0])/(1 + b*(y[3]*y[1] + (1 - y[3])*y[2])); //resource 1 
    dydt[3] = y[4]; //beta
    dydt[4] = v*p*y[0]*(y[1] - y[2])/pow(1 + b*(y[3]*y[1] + (1 - y[3])*y[2]) , 2); //gamma (helper function)

    return GSL_SUCCESS;
}

int function_chain(double t, const double y[], double dydt[], void* params){
    (void)(t);
    double* params_arr = (double*) params;
    double p = params_arr[0];
    double m = params_arr[1];
    double b = params_arr[2];
    double* alpha = &params_arr[3];
    double* beta = &params_arr[3 + N_RESOURCES*N_RESOURCES];
    dydt[0] = -m*y[0] + p*(y[3]*beta[0] + y[4]*beta[1] + y[5]*beta[2])*y[0]/(b*(y[3]*beta[0] + y[4]*beta[1] + y[5]*beta[2]) + 1);
    dydt[1] = -m*y[1] + p*(y[3]*beta[3] + y[4]*beta[4] + y[5]*beta[5])*y[1]/(b*(y[3]*beta[3] + y[4]*beta[4] + y[5]*beta[5]) + 1);
    dydt[2] = -m*y[2] + p*(y[3]*beta[6] + y[4]*beta[7] + y[5]*beta[8])*y[2]/(b*(y[3]*beta[6] + y[4]*beta[7] + y[5]*beta[8]) + 1);
    dydt[3] = (-p*y[2]*beta[6]/(b*(y[3]*beta[6] + y[4]*beta[7] + y[5]*beta[8]) + 1) - p*y[1]*beta[3]/(b*(y[3]*beta[3] + y[4]*beta[4] + y[5]*beta[5]) + 1) - p*y[0]*beta[0]/(b*(y[3]*beta[0] + y[4]*beta[1] + y[5]*beta[2]) + 1))*y[3] + (-y[3]*alpha[0] - y[4]*alpha[1] - y[5]*alpha[2] + 1)*y[3];
    dydt[4] = (-p*y[2]*beta[7]/(b*(y[3]*beta[6] + y[4]*beta[7] + y[5]*beta[8]) + 1) - p*y[1]*beta[4]/(b*(y[3]*beta[3] + y[4]*beta[4] + y[5]*beta[5]) + 1) - p*y[0]*beta[1]/(b*(y[3]*beta[0] + y[4]*beta[1] + y[5]*beta[2]) + 1))*y[4] + (-y[3]*alpha[3] - y[4]*alpha[4] - y[5]*alpha[5] + 1)*y[4];
    dydt[5] = (-p*y[2]*beta[8]/(b*(y[3]*beta[6] + y[4]*beta[7] + y[5]*beta[8]) + 1) - p*y[1]*beta[5]/(b*(y[3]*beta[3] + y[4]*beta[4] + y[5]*beta[5]) + 1) - p*y[0]*beta[2]/(b*(y[3]*beta[0] + y[4]*beta[1] + y[5]*beta[2]) + 1))*y[5] + (-y[3]*alpha[6] - y[4]*alpha[7] - y[5]*alpha[8] + 1)*y[5];

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

int jacobian_chain(double t, const double y[], double* dfdy, double dfdt[], void* params){
    (void)(t);
    double* params_arr = (double*) params;
    double p = params_arr[0];
    double m = params_arr[1];
    double b = params_arr[2];
    double* alpha = &params_arr[3];
    double* beta = &params_arr[3 + N_RESOURCES*N_RESOURCES];
    gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, N_EQUATIONS, N_EQUATIONS);
    gsl_matrix* jacobian_matrix = &dfdy_mat.matrix;
    gsl_matrix_set (jacobian_matrix, 0, 0, -m + p*(y[3]*beta[0] + y[4]*beta[1] + y[5]*beta[2])/(b*(y[3]*beta[0] + y[4]*beta[1] + y[5]*beta[2]) + 1));
    gsl_matrix_set (jacobian_matrix, 0, 1, 0);
    gsl_matrix_set (jacobian_matrix, 0, 2, 0);
    gsl_matrix_set (jacobian_matrix, 0, 3, -b*p*(y[3]*beta[0] + y[4]*beta[1] + y[5]*beta[2])*y[0]*beta[0]/pow(b*(y[3]*beta[0] + y[4]*beta[1] + y[5]*beta[2]) + 1, 2) + p*y[0]*beta[0]/(b*(y[3]*beta[0] + y[4]*beta[1] + y[5]*beta[2]) + 1));
    gsl_matrix_set (jacobian_matrix, 0, 4, -b*p*(y[3]*beta[0] + y[4]*beta[1] + y[5]*beta[2])*y[0]*beta[1]/pow(b*(y[3]*beta[0] + y[4]*beta[1] + y[5]*beta[2]) + 1, 2) + p*y[0]*beta[1]/(b*(y[3]*beta[0] + y[4]*beta[1] + y[5]*beta[2]) + 1));
    gsl_matrix_set (jacobian_matrix, 0, 5, -b*p*(y[3]*beta[0] + y[4]*beta[1] + y[5]*beta[2])*y[0]*beta[2]/pow(b*(y[3]*beta[0] + y[4]*beta[1] + y[5]*beta[2]) + 1, 2) + p*y[0]*beta[2]/(b*(y[3]*beta[0] + y[4]*beta[1] + y[5]*beta[2]) + 1));
    gsl_matrix_set (jacobian_matrix, 1, 0, 0);
    gsl_matrix_set (jacobian_matrix, 1, 1, -m + p*(y[3]*beta[3] + y[4]*beta[4] + y[5]*beta[5])/(b*(y[3]*beta[3] + y[4]*beta[4] + y[5]*beta[5]) + 1));
    gsl_matrix_set (jacobian_matrix, 1, 2, 0);
    gsl_matrix_set (jacobian_matrix, 1, 3, -b*p*(y[3]*beta[3] + y[4]*beta[4] + y[5]*beta[5])*y[1]*beta[3]/pow(b*(y[3]*beta[3] + y[4]*beta[4] + y[5]*beta[5]) + 1, 2) + p*y[1]*beta[3]/(b*(y[3]*beta[3] + y[4]*beta[4] + y[5]*beta[5]) + 1));
    gsl_matrix_set (jacobian_matrix, 1, 4, -b*p*(y[3]*beta[3] + y[4]*beta[4] + y[5]*beta[5])*y[1]*beta[4]/pow(b*(y[3]*beta[3] + y[4]*beta[4] + y[5]*beta[5]) + 1, 2) + p*y[1]*beta[4]/(b*(y[3]*beta[3] + y[4]*beta[4] + y[5]*beta[5]) + 1));
    gsl_matrix_set (jacobian_matrix, 1, 5, -b*p*(y[3]*beta[3] + y[4]*beta[4] + y[5]*beta[5])*y[1]*beta[5]/pow(b*(y[3]*beta[3] + y[4]*beta[4] + y[5]*beta[5]) + 1, 2) + p*y[1]*beta[5]/(b*(y[3]*beta[3] + y[4]*beta[4] + y[5]*beta[5]) + 1));
    gsl_matrix_set (jacobian_matrix, 2, 0, 0);
    gsl_matrix_set (jacobian_matrix, 2, 1, 0);
    gsl_matrix_set (jacobian_matrix, 2, 2, -m + p*(y[3]*beta[6] + y[4]*beta[7] + y[5]*beta[8])/(b*(y[3]*beta[6] + y[4]*beta[7] + y[5]*beta[8]) + 1));
    gsl_matrix_set (jacobian_matrix, 2, 3, -b*p*(y[3]*beta[6] + y[4]*beta[7] + y[5]*beta[8])*y[2]*beta[6]/pow(b*(y[3]*beta[6] + y[4]*beta[7] + y[5]*beta[8]) + 1, 2) + p*y[2]*beta[6]/(b*(y[3]*beta[6] + y[4]*beta[7] + y[5]*beta[8]) + 1));
    gsl_matrix_set (jacobian_matrix, 2, 4, -b*p*(y[3]*beta[6] + y[4]*beta[7] + y[5]*beta[8])*y[2]*beta[7]/pow(b*(y[3]*beta[6] + y[4]*beta[7] + y[5]*beta[8]) + 1, 2) + p*y[2]*beta[7]/(b*(y[3]*beta[6] + y[4]*beta[7] + y[5]*beta[8]) + 1));
    gsl_matrix_set (jacobian_matrix, 2, 5, -b*p*(y[3]*beta[6] + y[4]*beta[7] + y[5]*beta[8])*y[2]*beta[8]/pow(b*(y[3]*beta[6] + y[4]*beta[7] + y[5]*beta[8]) + 1, 2) + p*y[2]*beta[8]/(b*(y[3]*beta[6] + y[4]*beta[7] + y[5]*beta[8]) + 1));
    gsl_matrix_set (jacobian_matrix, 3, 0, -p*y[3]*beta[0]/(b*(y[3]*beta[0] + y[4]*beta[1] + y[5]*beta[2]) + 1));
    gsl_matrix_set (jacobian_matrix, 3, 1, -p*y[3]*beta[3]/(b*(y[3]*beta[3] + y[4]*beta[4] + y[5]*beta[5]) + 1));
    gsl_matrix_set (jacobian_matrix, 3, 2, -p*y[3]*beta[6]/(b*(y[3]*beta[6] + y[4]*beta[7] + y[5]*beta[8]) + 1));
    gsl_matrix_set (jacobian_matrix, 3, 3, -p*y[2]*beta[6]/(b*(y[3]*beta[6] + y[4]*beta[7] + y[5]*beta[8]) + 1) - p*y[1]*beta[3]/(b*(y[3]*beta[3] + y[4]*beta[4] + y[5]*beta[5]) + 1) - p*y[0]*beta[0]/(b*(y[3]*beta[0] + y[4]*beta[1] + y[5]*beta[2]) + 1) + (b*p*y[2]*pow(beta[6], 2)/pow(b*(y[3]*beta[6] + y[4]*beta[7] + y[5]*beta[8]) + 1, 2) + b*p*y[1]*pow(beta[3], 2)/pow(b*(y[3]*beta[3] + y[4]*beta[4] + y[5]*beta[5]) + 1, 2) + b*p*y[0]*pow(beta[0], 2)/pow(b*(y[3]*beta[0] + y[4]*beta[1] + y[5]*beta[2]) + 1, 2))*y[3] - 2*y[3]*alpha[0] - y[4]*alpha[1] - y[5]*alpha[2] + 1);
    gsl_matrix_set (jacobian_matrix, 3, 4, (b*p*y[2]*beta[6]*beta[7]/pow(b*(y[3]*beta[6] + y[4]*beta[7] + y[5]*beta[8]) + 1, 2) + b*p*y[1]*beta[3]*beta[4]/pow(b*(y[3]*beta[3] + y[4]*beta[4] + y[5]*beta[5]) + 1, 2) + b*p*y[0]*beta[0]*beta[1]/pow(b*(y[3]*beta[0] + y[4]*beta[1] + y[5]*beta[2]) + 1, 2))*y[3] - y[3]*alpha[1]);
    gsl_matrix_set (jacobian_matrix, 3, 5, (b*p*y[2]*beta[6]*beta[8]/pow(b*(y[3]*beta[6] + y[4]*beta[7] + y[5]*beta[8]) + 1, 2) + b*p*y[1]*beta[3]*beta[5]/pow(b*(y[3]*beta[3] + y[4]*beta[4] + y[5]*beta[5]) + 1, 2) + b*p*y[0]*beta[0]*beta[2]/pow(b*(y[3]*beta[0] + y[4]*beta[1] + y[5]*beta[2]) + 1, 2))*y[3] - y[3]*alpha[2]);
    gsl_matrix_set (jacobian_matrix, 4, 0, -p*y[4]*beta[1]/(b*(y[3]*beta[0] + y[4]*beta[1] + y[5]*beta[2]) + 1));
    gsl_matrix_set (jacobian_matrix, 4, 1, -p*y[4]*beta[4]/(b*(y[3]*beta[3] + y[4]*beta[4] + y[5]*beta[5]) + 1));
    gsl_matrix_set (jacobian_matrix, 4, 2, -p*y[4]*beta[7]/(b*(y[3]*beta[6] + y[4]*beta[7] + y[5]*beta[8]) + 1));
    gsl_matrix_set (jacobian_matrix, 4, 3, (b*p*y[2]*beta[6]*beta[7]/pow(b*(y[3]*beta[6] + y[4]*beta[7] + y[5]*beta[8]) + 1, 2) + b*p*y[1]*beta[3]*beta[4]/pow(b*(y[3]*beta[3] + y[4]*beta[4] + y[5]*beta[5]) + 1, 2) + b*p*y[0]*beta[0]*beta[1]/pow(b*(y[3]*beta[0] + y[4]*beta[1] + y[5]*beta[2]) + 1, 2))*y[4] - y[4]*alpha[3]);
    gsl_matrix_set (jacobian_matrix, 4, 4, -p*y[2]*beta[7]/(b*(y[3]*beta[6] + y[4]*beta[7] + y[5]*beta[8]) + 1) - p*y[1]*beta[4]/(b*(y[3]*beta[3] + y[4]*beta[4] + y[5]*beta[5]) + 1) - p*y[0]*beta[1]/(b*(y[3]*beta[0] + y[4]*beta[1] + y[5]*beta[2]) + 1) + (b*p*y[2]*pow(beta[7], 2)/pow(b*(y[3]*beta[6] + y[4]*beta[7] + y[5]*beta[8]) + 1, 2) + b*p*y[1]*pow(beta[4], 2)/pow(b*(y[3]*beta[3] + y[4]*beta[4] + y[5]*beta[5]) + 1, 2) + b*p*y[0]*pow(beta[1], 2)/pow(b*(y[3]*beta[0] + y[4]*beta[1] + y[5]*beta[2]) + 1, 2))*y[4] - y[3]*alpha[3] - 2*y[4]*alpha[4] - y[5]*alpha[5] + 1);
    gsl_matrix_set (jacobian_matrix, 4, 5, (b*p*y[2]*beta[7]*beta[8]/pow(b*(y[3]*beta[6] + y[4]*beta[7] + y[5]*beta[8]) + 1, 2) + b*p*y[1]*beta[4]*beta[5]/pow(b*(y[3]*beta[3] + y[4]*beta[4] + y[5]*beta[5]) + 1, 2) + b*p*y[0]*beta[1]*beta[2]/pow(b*(y[3]*beta[0] + y[4]*beta[1] + y[5]*beta[2]) + 1, 2))*y[4] - y[4]*alpha[5]);
    gsl_matrix_set (jacobian_matrix, 5, 0, -p*y[5]*beta[2]/(b*(y[3]*beta[0] + y[4]*beta[1] + y[5]*beta[2]) + 1));
    gsl_matrix_set (jacobian_matrix, 5, 1, -p*y[5]*beta[5]/(b*(y[3]*beta[3] + y[4]*beta[4] + y[5]*beta[5]) + 1));
    gsl_matrix_set (jacobian_matrix, 5, 2, -p*y[5]*beta[8]/(b*(y[3]*beta[6] + y[4]*beta[7] + y[5]*beta[8]) + 1));
    gsl_matrix_set (jacobian_matrix, 5, 3, (b*p*y[2]*beta[6]*beta[8]/pow(b*(y[3]*beta[6] + y[4]*beta[7] + y[5]*beta[8]) + 1, 2) + b*p*y[1]*beta[3]*beta[5]/pow(b*(y[3]*beta[3] + y[4]*beta[4] + y[5]*beta[5]) + 1, 2) + b*p*y[0]*beta[0]*beta[2]/pow(b*(y[3]*beta[0] + y[4]*beta[1] + y[5]*beta[2]) + 1, 2))*y[5] - y[5]*alpha[6]);
    gsl_matrix_set (jacobian_matrix, 5, 4, (b*p*y[2]*beta[7]*beta[8]/pow(b*(y[3]*beta[6] + y[4]*beta[7] + y[5]*beta[8]) + 1, 2) + b*p*y[1]*beta[4]*beta[5]/pow(b*(y[3]*beta[3] + y[4]*beta[4] + y[5]*beta[5]) + 1, 2) + b*p*y[0]*beta[1]*beta[2]/pow(b*(y[3]*beta[0] + y[4]*beta[1] + y[5]*beta[2]) + 1, 2))*y[5] - y[5]*alpha[7]);
    gsl_matrix_set (jacobian_matrix, 5, 5, -p*y[2]*beta[8]/(b*(y[3]*beta[6] + y[4]*beta[7] + y[5]*beta[8]) + 1) - p*y[1]*beta[5]/(b*(y[3]*beta[3] + y[4]*beta[4] + y[5]*beta[5]) + 1) - p*y[0]*beta[2]/(b*(y[3]*beta[0] + y[4]*beta[1] + y[5]*beta[2]) + 1) + (b*p*y[2]*pow(beta[8], 2)/pow(b*(y[3]*beta[6] + y[4]*beta[7] + y[5]*beta[8]) + 1, 2) + b*p*y[1]*pow(beta[5], 2)/pow(b*(y[3]*beta[3] + y[4]*beta[4] + y[5]*beta[5]) + 1, 2) + b*p*y[0]*pow(beta[2], 2)/pow(b*(y[3]*beta[0] + y[4]*beta[1] + y[5]*beta[2]) + 1, 2))*y[5] - y[3]*alpha[6] - y[4]*alpha[7] - 2*y[5]*alpha[8] + 1);

    for (int i = 0; i < N_EQUATIONS; i++){
        dfdt[i] = 0.;
    }
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

bool check_beta_convergence(const double& beta_conv, const double* y, double y_prev[][2], double dt){

    if (abs(y[3]) < beta_conv || abs(y[3] - 1) < beta_conv ){
        beta_convergence_counter += dt;
    }
    else {
        beta_convergence_counter = 0;   
    }

    if (beta_convergence_counter >= BETA_CONVERGENCE_TIME){
        return true;
    }
    return false;
}

void check_extremum(const double* y, double y_prev[][2], const double& t){
    for (int i = 0; i < N_EQUATIONS; i++){
        if (y_prev[i][0] > y[i] && y_prev[i][0] > y_prev[i][1]){
            maxima.at(i).at(maxima_counter.at(i)%N_EXTREMA) = y_prev[i][0];
            maxima_time.at(i).at(maxima_counter.at(i)%N_EXTREMA) = t;
            //std::cout << "Maxima at " << i << "\t" << maxima.at(i).at(maxima_counter.at(i)) << std::endl; 
            maxima_counter.at(i) ++;
        }
        else if(y_prev[i][0] < y[i] && y_prev[i][0] < y_prev[i][1]){
            minima.at(i).at(minima_counter.at(i)%N_EXTREMA) = y_prev[i][0];
            minima_time.at(i).at(minima_counter.at(i)%N_EXTREMA) = t;
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

void save_periods(std::ofstream& OutputPeriods){
    for (int i = 0 ; i<N_EQUATIONS; i++){
        double period {};
        for (int j = 0 ; j < N_EXTREMA - 1; j++){
            period += abs(maxima_time.at(i).at(j+1) - maxima_time.at(i).at(j));
            period += abs(minima_time.at(i).at(j+1) - minima_time.at(i).at(j));
        }
        OutputPeriods << period/(2*N_EXTREMA) << std::endl;
    }
}

void save_node(std::ofstream& OutputNode, double* y){
    for (int i = 0; i < N_EQUATIONS; i++){
        OutputNode << y[i];
        OutputNode << "\t"; //gurantees NULL at reading second column
        OutputNode << std::endl;
    }
}
void save_time_dependence(std::ofstream& OutputTimeDependence, const double* y, const double& t){
    OutputTimeDependence << t;
    for (int i = 0; i < N_EQUATIONS; i++){
        OutputTimeDependence << '\t' << y[i];
    }
    OutputTimeDependence << std::endl;
}

void solve_system(double dt, double t_max, double t0, double* y0, double eps_abs, double eps_rel, void* params, \
                const double& conv, std::ofstream& OutputTimeDependence, std::ofstream& OutputCycleOrNode, std::ofstream& OutputPeriods, f_equation_system system, f_jacobian jacobian){
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
        //OutputTimeDependence << t << "\t" << y[0] << "\t" << y[1] << "\t" << y[2] << "\t" << y[3] << "\t" << y[4] << std::endl;
        save_time_dependence(OutputTimeDependence, y, t);

        if (status != GSL_SUCCESS) {
            printf ("error: driver returned %d at time %f\n", status, t);
            break;
        }

        //Check if beta [0,1]
        // if (y[3] < 0)
        //     y[3] = 0.;
        // else if (y[3] > 1)
        //     y[3] = 1;

        check_extremum(y, y_prev, t);

        conv_flag = check_convergence(conv, y, y_prev);
        //conv_flag = check_beta_convergence(conv, y, y_prev, dt);
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
        save_periods(OutputPeriods);
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
    double eps_abs = 1e-10;
    double eps_rel = 1e-10;
    double convergence = 1e-8;

    double y0[N_EQUATIONS] = {};
    // params = {p , m , b, alpha_11, alpha_12, alpha_21, alpha_22, v}
    double params[N_PARAMS] {};

    double dt = std::stod(argv[1]);
    double t_max = std::stod(argv[2]);
    for (int i = 0; i < N_EQUATIONS; i++){
        y0[i] = std::stod(argv[3 + i]);
    }

    for (int i = 0; i < N_PARAMS; i++){
        params[i] = std::stod(argv[3 + N_EQUATIONS + i]);
    }
    std::string path = argv[3 + N_EQUATIONS + N_PARAMS];


    std::string filename = path + "/TimeDependence.dat";
    std::ofstream OutputTimeDependence(filename);
    std::string filename_final = path + "/FinalState.dat";
    std::ofstream OutputCycleOrNode(filename_final);
    std::string filename_period = path + "/Periods.dat";
    std::ofstream OutputPeriods(filename_period);
    solve_system(dt, t_max, t0, y0, eps_abs, eps_rel, params, convergence, OutputTimeDependence, OutputCycleOrNode, OutputPeriods, function_chain, jacobian_chain);

    return 0;
}
