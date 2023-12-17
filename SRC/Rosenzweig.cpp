#include "Equations.h"

inline double rand_uniform(void){ return double(rand())/RAND_MAX;}

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
            maxima_time.at(i) += t - prev_maximum.at(i);
            prev_maximum.at(i) = t;
            //std::cout << "Equation: " << i << "time: " << t << std::endl;
            //std::cout << "Maxima at " << i << "\t" << maxima.at(i).at(maxima_counter.at(i)) << std::endl; 
            maxima_counter.at(i) ++;
        }
        else if(y_prev[i][0] < y[i] && y_prev[i][0] < y_prev[i][1]){
            minima.at(i).at(minima_counter.at(i)%N_EXTREMA) = y_prev[i][0];
            minima_time.at(i) += t - prev_minimum.at(i);
            prev_minimum.at(i) = t;
            //std::cout << "Minima at " << i << "\t" << minima.at(i).at(minima_counter.at(i)) << std::endl; 
            minima_counter.at(i) ++;
        }
    }
    return;
}

void update_beta(void* params){
    double* params_arr = (double*) params;
    double* beta = &params_arr[3 + N_RESOURCES*N_RESOURCES];
    size_t beta_size = N_RESOURCES*N_CONSUMERS;
    double max_beta {};
    for (size_t i = 0; i < beta_size; i++){
        if (beta[i] > max_beta) max_beta = beta[i];
    }

    for (size_t i = 0; i < N_CONSUMERS; i++){
        if (rand_uniform() > 0.5){
            beta[i*N_RESOURCES + i] = max_beta;
            beta[i*N_RESOURCES + (i + 1)%N_RESOURCES] = max_beta;
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
        period = (maxima_time.at(i) / maxima_counter.at(i) + minima_time.at(i) / minima_counter.at(i))/2.;
        OutputPeriods << period << std::endl;
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
                const double& conv, const double& T_beta_update, std::ofstream& OutputTimeDependence, std::ofstream& OutputCycleOrNode, std::ofstream& OutputPeriods, f_equation_system system, f_jacobian jacobian){
    //GSL procedures allocation
    gsl_odeiv2_system my_equation = {system, jacobian, N_EQUATIONS, params};
    gsl_odeiv2_driver* driver = gsl_odeiv2_driver_alloc_y_new (&my_equation, gsl_odeiv2_step_rk8pd, dt, eps_abs, eps_rel);
    
    size_t n_updates = 1;

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

        if (t > T_beta_update*n_updates){
            update_beta(params);
            n_updates ++;
        }


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
        if (argc != 5 + N_EQUATIONS + N_PARAMS){
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
    double T_beta_update = std::stod(argv[3]);
    for (int i = 0; i < N_EQUATIONS; i++){
        y0[i] = std::stod(argv[4 + i]);
    }

    for (int i = 0; i < N_PARAMS; i++){
        params[i] = std::stod(argv[4 + N_EQUATIONS + i]);
    }
    std::string path = argv[4 + N_EQUATIONS + N_PARAMS];


    std::string filename = path + "/TimeDependence.dat";
    std::ofstream OutputTimeDependence(filename);
    std::string filename_final = path + "/FinalState.dat";
    std::ofstream OutputCycleOrNode(filename_final);
    std::string filename_period = path + "/Periods.dat";
    std::ofstream OutputPeriods(filename_period);
    solve_system(dt, t_max, t0, y0, eps_abs, eps_rel, params, convergence, T_beta_update, OutputTimeDependence, OutputCycleOrNode, OutputPeriods, function_chain, jacobian_chain);

    return 0;
}
