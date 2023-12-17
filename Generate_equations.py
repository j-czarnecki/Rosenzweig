import numpy as np
from sympy import *
from sympy.utilities.codegen import codegen
import os


###################################### SYMBOLIC OPERATIONS #################################
def generate_equations(n_consumers, n_resources, output_path):
    print("Generating equations for n_consumers = " + str(n_consumers) + " n_resources = " + str(n_resources))
    #n_consumers = 3
    #n_resources = 3
    t = Symbol('t')

    #first we put n_consumers to y, then n_resources
    y_symbols = symbols('y:{0}'.format(n_consumers + n_resources), cls = Function)
    y = [yi(t) for yi in y_symbols]
    #print(y)

    #matrix defining interactions between resources
    alpha_mat = Matrix(MatrixSymbol('alpha', n_resources, n_resources))
    #print(alpha_mat)
    #matrix defining interactions between consumers and resources
    #i.e. beta_mat[i,j] is coefficient of i-th consumer taking advantage of j-th  resource
    beta_mat = Matrix(MatrixSymbol('beta', n_consumers, n_resources))
    #print(beta_mat)
    #m - mortality rate, p - hunting rate and b - Holling constant
    m, p, b = symbols('m p b')
    #print(m,p,b)


    omega = 0
    for i in range(n_consumers):
        for j in range(n_resources):
            omega += beta_mat[i,j]*y[n_consumers + j]

    omega = omega*b + 1
    #print(omega)

    f = []
    for i in range(n_consumers + n_resources):

        fi = 0
        if(i < n_consumers): #consumers equations

            hunting_impact = 0
            for j in range(n_resources):
                hunting_impact += beta_mat[i, j]*y[n_consumers + j]
            hunting_impact *= p*y[i]/(1 + b*hunting_impact)     
            
            fi += hunting_impact - m*y[i]

        elif(i >= n_consumers): #resources equations
            other_resources_impact = 1
            for j in range(n_resources):
                other_resources_impact -= alpha_mat[i - n_consumers,j]*y[n_consumers + j] 
            other_resources_impact = y[i]*other_resources_impact

            consumers_impact = 0
            for j in range(n_consumers):
                single_consumer_impact = p*y[j]*beta_mat[j, i - n_consumers]
                hunting_impact = 0
                for k in range(n_resources):
                    hunting_impact += beta_mat[j,k]*y[n_consumers + k]            
                consumers_impact -= single_consumer_impact/(1 + b*hunting_impact)
            consumers_impact *= y[i]     
            
            fi = other_resources_impact + consumers_impact

        f.append(fi)
    #print('\n\n')
    #for i in range(n_consumers + n_resources):
        #print(f[i])

    #print('\n\n')
    jacobian = Matrix([[diff(fi, variable) for variable in y] for fi in f])
    ###################### END OF SYMBOLIC OPERATIONS #################################################


    #################### GSL CODE GENERATION ##########################################################

    ccode_jacobian = ""
    #print jacobian
    for i in range(n_consumers + n_resources):
        for j in range(n_consumers + n_resources):
            elem = jacobian[i,j]
            #returns a tuple with elements (symbols, errors, code)
            #errors are dealt with during parsing
            ccode_elem = ccode(elem, human = False)[2] 
            gsl_matrix_assignment = f"\tgsl_matrix_set (jacobian_matrix, {i}, {j}, ".format(i,j)
            ccode_elem_gsl = gsl_matrix_assignment + ccode_elem + ");\n"
            ccode_jacobian += ccode_elem_gsl

    #change y0(t) -> y[0] etc.
    for i in range(n_consumers + n_resources):
        current_y = f"y{i}(t)".format(i)
        ccode_y = f"y[{i}]".format(i)
        ccode_jacobian = ccode_jacobian.replace(current_y, ccode_y)
    #print jacobian for gsl to file 
    with open(os.path.join(output_path, "Jacobian.ccode"), "w") as jacobian_file:
        print(ccode_jacobian, file = jacobian_file)



    #print('\n\n')
    #print equations
    ccode_equation_system = ""
    for i in range(n_consumers + n_resources):
        ccode_equation_system += '\t' + ccode(f[i], assign_to=f"dydt[{i}]".format(i), standard='C99', human = False)[2] + '\n'

    #change y0(t) -> y[0] etc.
    for i in range(n_consumers + n_resources):
        current_y = f"y{i}(t)".format(i)
        ccode_y = f"y[{i}]".format(i)
        ccode_equation_system = ccode_equation_system.replace(current_y, ccode_y)
    #print(ccode_equation_system)

    ###################################### Equations_constants.h ###########################################
    equations_constants_header = "#include <iostream> \n#include <cmath>\n#include <gsl/gsl_errno.h>\n#include <gsl/gsl_matrix.h> \
    \n#include <gsl/gsl_odeiv2.h>\n#include <fstream>\n#include <vector>\n#include <array>\n#include <numeric>\n"

    global_data = "//Global data initialization\n\
    const int N_CONSUMERS = " + str(n_consumers) + ";\n\
    const int N_RESOURCES = " + str(n_resources) + ";\n\
    const int N_EQUATIONS = N_CONSUMERS + N_RESOURCES;\n\
    const int N_PARAMS = 3 + N_RESOURCES*N_RESOURCES + N_CONSUMERS*N_RESOURCES;\n\
    const int N_EXTREMA = 100;\n\
    const double BETA_CONVERGENCE_TIME= 5e4;\n\
    \n\
    double beta_convergence_counter = 0;\n\
    // std::vector<std::vector<double>> maxima(N_EQUATIONS);\n\
    // std::vector<std::vector<double>> minima(N_EQUATIONS);\n\
    std::array<std::array<double, N_EXTREMA>, N_EQUATIONS> maxima;\n\
    std::array<std::array<double, N_EXTREMA>, N_EQUATIONS> minima;\n\
    std::array<double, N_EQUATIONS> maxima_time;\n\
    std::array<double, N_EQUATIONS> minima_time;\n\
    std::array<double, N_EQUATIONS> prev_maximum;\n\
    std::array<double, N_EQUATIONS> prev_minimum;\n\
    std::array<size_t, N_EQUATIONS> maxima_counter {};\n\
    std::array<size_t, N_EQUATIONS> minima_counter {};\n"

    with open(os.path.join(output_path, "Equations_constants.h"), "w") as constants_file:
        print(equations_constants_header, file = constants_file)
        print(global_data, file = constants_file)


    ########################################## Equations.h ########################################
    include_constants_header = "#include \"Equations_constants.h\"\n"

    typedefs = "typedef int (*f_equation_system)(double, const double*, double*, void*);\n\
    typedef int (*f_jacobian)(double, const double*, double*, double*, void*);\n"

    #function
    function_declaration = "int function_chain(double t, const double y[], double dydt[], void* params){"
    unpack_params = "\t(void)(t);\n\
    \tdouble* params_arr = (double*) params;\n\
    \tdouble p = params_arr[0];\n\
    \tdouble m = params_arr[1];\n\
    \tdouble b = params_arr[2];\n\
    \tdouble* alpha = &params_arr[3];\n\
    \tdouble* beta = &params_arr[3 + N_RESOURCES*N_RESOURCES];\n"
    return_gsl = "  return GSL_SUCCESS;\n"

    #jacobian
    jacobian_declaration = "int jacobian_chain(double t, const double y[], double* dfdy, double dfdt[], void* params){"
    gsl_matrix = "\tgsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, N_EQUATIONS, N_EQUATIONS);\n\
    \tgsl_matrix* jacobian_matrix = &dfdy_mat.matrix;\n"

    time_partial = "\tfor (int i = 0; i < N_EQUATIONS; i++){\n\
    \t\tdfdt[i] = 0.;\n\
    \t}"


    print("Saving functions to header file")
    with open(os.path.join(output_path, "Equations.h"), "w") as equations_file:
        print(include_constants_header, file = equations_file)
        print(typedefs, file=equations_file)
        #Equations
        print(function_declaration, file = equations_file)
        print(unpack_params, file = equations_file)
        print(ccode_equation_system, file = equations_file)
        print(return_gsl, file = equations_file)
        print("}\n", file = equations_file)
        #Jacobian
        print(jacobian_declaration, file = equations_file)
        print(unpack_params, file = equations_file)
        print(gsl_matrix, file = equations_file)
        print(ccode_jacobian, file = equations_file)
        print(time_partial, file = equations_file)
        print(return_gsl, file = equations_file)
        print("}\n", file = equations_file)

    print("Header created")

