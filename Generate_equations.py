import numpy as np
from sympy import *
from sympy.utilities.codegen import codegen


###################################### SYMBOLIC OPERATIONS #################################
n_consumers = 3
n_resources = 3
t = Symbol('t')

#first we put n_consumers to y, then n_resources
y_symbols = symbols('y:{0}'.format(n_consumers + n_resources), cls = Function)
y = [yi(t) for yi in y_symbols]
print(y)

#matrix defining interactions between resources
alpha_mat = Matrix(MatrixSymbol('alpha', n_resources, n_resources))
print(alpha_mat)
#matrix defining interactions between consumers and resources
#i.e. beta_mat[i,j] is coefficient of i-th consumer taking advantage of j-th  resource
beta_mat = Matrix(MatrixSymbol('beta', n_consumers, n_resources))
print(beta_mat)
#m - mortality rate, p - hunting rate and b - Holling constant
m, p, b = symbols('m p b')
print(m,p,b)


omega = 0
for i in range(n_consumers):
    for j in range(n_resources):
        omega += beta_mat[i,j]*y[n_consumers + j]

omega = omega*b + 1
print(omega)

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
print('\n\n')
for i in range(n_consumers + n_resources):
    print(f[i])

print('\n\n')
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
        gsl_matrix_assignment = f"gsl_matrix_set (jacobian_matrix, {i}, {j}, ".format(i,j)
        ccode_elem_gsl = gsl_matrix_assignment + ccode_elem + ");\n"
        ccode_jacobian += ccode_elem_gsl

#change y0(t) -> y[0] etc.
for i in range(n_consumers + n_resources):
    current_y = f"y{i}(t)".format(i)
    ccode_y = f"y[{i}]".format(i)
    ccode_jacobian = ccode_jacobian.replace(current_y, ccode_y)
#print jacobian for gsl to file 
with open("./Jacobian.ccode", "w") as jacobian_file:
    print(ccode_jacobian, file = jacobian_file)



print('\n\n')
#print equations
ccode_equation_system = ""
for i in range(n_consumers + n_resources):
    ccode_equation_system += ccode(f[i], assign_to=f"dydt[{i}]".format(i), standard='C99', human = False)[2] + '\n'

#change y0(t) -> y[0] etc.
for i in range(n_consumers + n_resources):
    current_y = f"y{i}(t)".format(i)
    ccode_y = f"y[{i}]".format(i)
    ccode_equation_system = ccode_equation_system.replace(current_y, ccode_y)
print(ccode_equation_system)
with open("./Equations.ccode", "w") as equations_file:
    print(ccode_equation_system, file = equations_file)


