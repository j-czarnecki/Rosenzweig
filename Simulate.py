import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import pandas as pd
import numpy as np
import os 
import subprocess
import matplotlib.collections as mcoll
import matplotlib.path as mpath
import math
import random 
import Generate_equations
import multiprocessing
import time


N_EQUATIONS = 5
cwd = os.getcwd()
# alpha_tab = [0.5, 0.6, 0.7, 0.8, 0.9, 1.]
# beta_tab = [0.2, 0.4, 0.6, 0.8, 0.9, 1.]
#mortality_tab = [0.05*i for i in range(1,21)]

alpha_tab = [1., 0.98, 0.95, 0.9, 0.85, 0.8]
#alpha_tab = [.95]
v_tab = [10e-8*i for i in range(1, 41)]
final_states = [[0]*len(v_tab) for i in range(2*N_EQUATIONS)]
final_phase = [0 for i in range(len(v_tab))]

fig, ax = plt.subplots()


def clear_final_states():
    final_states = [[0]*len(v_tab) for i in range(2*N_EQUATIONS)]

def clear_final_phase():
    final_phase = [0 for i in range(len(v_tab))]


def plot_phase(run_cmd):
    path = run_cmd[16]
    alpha_11 = float(run_cmd[11])
    v = float(run_cmd[15])
    column_names = ["t", "C", "R1", "R2", "beta", "gamma"]

    data = pd.read_table(f"{path}/TimeDependence.dat", delimiter="\t", names=column_names)    
    plt.clf()
    plt.plot(data.R1, data.C)
    plt.xlabel(r"$R_1$")
    plt.ylabel(r"C")
    plt.grid(1)
    plt.savefig(f"{cwd}/Plots/PhaseCR1_alpha{alpha_11:.3f}_v{v*1e8:.3f}.png")

    plt.clf()
    plt.plot(data.R2, data.C)
    plt.xlabel(r"$R_2$")
    plt.ylabel(r"C")
    plt.grid(1)
    plt.savefig(f"{cwd}/Plots/PhaseCR2_alpha{alpha_11:.3f}_v{v*1e8:.3f}.png")

    plt.clf()
    plt.plot(data.R1, data.R2)
    plt.xlabel(r"$R_1$")
    plt.ylabel(r"$R_2$")
    plt.grid(1)
    plt.savefig(f"{cwd}/Plots/PhaseR2R1_alpha{alpha_11:.3f}_v{v*1e8:.3f}.png")

    plt.clf()
    plt.plot(data.gamma, data.beta)
    plt.xlabel(r"$\gamma$")
    plt.ylabel(r"$\beta$")
    plt.grid(1)
    plt.savefig(f"{cwd}/Plots/PhaseBetaGamma_alpha{alpha_11:.3f}_v{v*1e8:.3f}.png")


def plot_time_dependence(run_cmd):
    path = run_cmd[16]
    alpha_11 = float(run_cmd[11])
    v = float(run_cmd[15])

    column_names = ["t", "C", "R1", "R2", "beta", "gamma"]
    labels = ["t", "C", r"$R_1$", r"$R_2$", r"$ \beta $", r"$\gamma$"]
    data = pd.read_table(f"{path}/TimeDependence.dat", delimiter="\t", names=column_names)
    for i in range(1, len(column_names)):
        plt.clf()
        plt.plot(data.t, data[column_names[i]])
        plt.xlabel("t")
        plt.ylabel(labels[i])
        plt.grid(1)
        plt.tight_layout()
        plt.savefig(f"{cwd}/Plots/{column_names[i]}(t)_alpha{alpha_11:.3f}_v{v*1e8:.3f}.png")


def plot_final_states(run_cmd):
    alpha_11 = float(run_cmd[11])

    styles = ['r.', 'k.', 'c.', 'g.', 'y.']
    labels = ["C", None, r"$R_1$", None, r"$R_2$", None, r"$ \beta $",None, r"$\gamma$", None]

    plt.clf()
    for i in range(2*3):
        plt.plot(v_tab, final_states[i][:], styles[int(i/2)], label = labels[i])
    plt.legend()
    plt.xlabel("v")
    #plt.ylabel("R,C")
    plt.grid(1)
    plt.savefig(f"{cwd}/Plots/FinalState_alpha{alpha_11:.3f}.png")

def plot_final_phase(run_cmd):
    alpha_11 = float(run_cmd[11])
    plt.clf()
    plt.plot(v_tab, final_phase)
    plt.xlabel("v")
    #plt.ylabel("R,C")
    plt.grid(1)
    plt.savefig(f"{cwd}/Plots/FinalPhase_alpha{alpha_11:.3f}.png")



def get_final_states(run_cmd, m_pos):
    path = run_cmd[16]
    data = pd.read_table(f"{path}/FinalState.dat", delimiter="\t", names=["maximum", "minimum"])
    for i in range(len(data.maximum)):
        final_states[2*i][m_pos] = data.maximum[i]
        final_states[2*i+1][m_pos] = data.minimum[i]

def get_final_phase(run_cmd, m_pos):
    path = run_cmd[16]
    data = pd.read_table(f"{path}/FinalState.dat", delimiter="\t", names=["maximum", "minimum"])
    print(data.minimum[0])
    if math.isnan(data.minimum[0]):
        final_phase[m_pos] = 'Stable'
        #print('Phase G found')
    else: 
        final_phase[m_pos] = 'Unstable'


def single_run(n_cells, output_path):
    print(f"Inside single run")
    n_consumers = n_cells
    n_resources = n_cells

    base_directory = os.getcwd()
    if not os.path.exists(output_path):
        os.mkdir(output_path)  

    SRC_dir = os.path.join(base_directory, "SRC")
    cp_arguments = ["cp", "-r", SRC_dir, os.path.join(base_directory,output_path)]
    cp_src = subprocess.run(cp_arguments)
    #make sure that previous process managed to run
    #with proper Equations.h file
    Generate_equations.generate_equations(n_consumers, n_resources, os.path.join(output_path, "SRC"))
    make_cmd = ["make", "-C", os.path.join(output_path, "SRC")]
    make = subprocess.run(make_cmd) #Always after generating equations new .x should be built

    dt = 1e-1
    t_max = 1e6
    p = 3.
    m = 0.2
    b = 2.
    T_beta_update = 1e10
    #output_path = f"OutputData_N" + str(n_cells)

    C0 = [0.4]*n_consumers #all consumers start with same population
    R0 = []
    for i in range(n_resources):
        R0.append(0.7 + random.uniform(-0.05, 0.05)) # Slightly varying resource population
    
    #Interaction between resources
    # alpha_ij denotes interaction between i-th and j-th resource
    alpha = np.zeros((n_resources, n_resources))
    for i in range(n_resources):
        alpha[i,i] = 0.95
        alpha[i,(i+1)%n_resources] = 0.065
        alpha[(i+1)%n_resources, i] = 0.065
    #To delete interacion between boundary resources
    #We have to do this, since with periodic BC
    #This is the same resource
        
    alpha[n_resources - 1,0] = alpha[0, n_resources - 1] = 0
    

    #Interaction between consumers nad resources
    #Beta matrix defines coefficient of i-th consumer taking advantage of j-th resource
    beta = np.zeros((n_consumers, n_resources))
    for i in range(n_consumers):
        beta[i,i] = np.random.uniform(0.475,0.525)
        beta[i, (i+1)%n_resources] = beta[i,i]
    beta[n_consumers - 1, 0] = 0.

    with open(os.path.join(output_path, "beta.dat"), "w") as beta_file:
        print("#Beta matrix defines coefficient of i-th consumer taking advantage of j-th resource", file = beta_file)
        print(beta, file = beta_file)
    
    with open(os.path.join(output_path, "alpha.dat"), "w") as alpha_file:
        print("#Alpha matrix defines interactions between i-th and j-th resource", file = alpha_file)
        print(alpha, file = alpha_file)

    #This is needed to easily append to argument list
    alpha = alpha.flatten()
    beta = beta.flatten()

    #Creating arguments list.
    #Order is crucial, since it is passed via argv to C++ program
    arguments = []
    arguments.append(str(dt))
    arguments.append(str(t_max))
    arguments.append(str(T_beta_update))
    for i in range(n_consumers):
        arguments.append(str(C0[i]))
    for i in range(n_resources):
        arguments.append(str(R0[i]))
    arguments.append(str(p))
    arguments.append(str(m))
    arguments.append(str(b))
    for i in range(len(alpha)):
        arguments.append(str(alpha[i]))
    for i in range(len(beta)):
        arguments.append(str(beta[i]))
    arguments.append(output_path)

    if not os.path.exists(arguments[-1]):
        os.mkdir(arguments[-1]) 
    run_with_arguments = [os.path.join(output_path, "SRC", 'Rosenzweig.x')] + arguments
    print(run_with_arguments)
    simulation = subprocess.run(run_with_arguments)

    return simulation


if __name__ == "__main__":
    #First element specifies number of chain cells (Resource - Consumer pairs)
    #Second element of tuple specifies output directory
    #N processes specifies how many processes can run simultaneously    
    input_tab = [(i, f"/mnt/HDD1TB/Rosenzweig/RUN_n" + str(i) + "_Stiff") for i in [3,4,5,10,25]]
    N_processes = 5

    with multiprocessing.Pool(processes=N_processes) as pool:
        result = pool.starmap(single_run, input_tab)
# def main():
#     dt = 1e-1
#     t_max = 2e6
#     C0 = 0.4
#     R1_0 = 0.71
#     R2_0 = 0.7 
#     beta0 = 0.49
#     gamma0 = 0.
#     p = 3.
#     m = 0.2
#     b = 2.
#     alpha_11 = 0.95
#     alpha_12 = 0.065
#     alpha_21 = 0.065
#     alpha_22 = 0.95
#     v = 80e-8

#     for alpha_11 in alpha_tab:
#         alpha_22 = alpha_11
#         v_pos = 0
#         for v in v_tab:
#             current_workdir = os.getcwd()
#             path = f"/mnt/HDD1TB/Rosenzweig/RUN_alpha{alpha_11:.3f}_v{v*1e8:.3f}"
#             if not os.path.exists(path):
#                 os.mkdir(path) 

#             run_with_arguments = [current_workdir + '/a.out'] + [str(dt), str(t_max), str(C0),\
#                                             str(R1_0), str(R2_0), str(beta0), str(gamma0), str(p), str(m), str(b), \
#                                             str(alpha_11), str(alpha_12), str(alpha_21), str(alpha_22), str(v), path]
#             print(run_with_arguments)
#             simulation = subprocess.run(run_with_arguments)
#             plot_time_dependence(run_with_arguments)
#             plot_phase(run_with_arguments)
#             #get_final_states(run_with_arguments, v_pos)
#             get_final_phase(run_with_arguments, v_pos)
#             v_pos += 1
#         #plot_final_states(run_with_arguments)
#         #clear_final_states()
#         plot_final_phase(run_with_arguments)
#         clear_final_phase()






# def plot_phase_C1R1(run_cmd):
#     alpha_11 = float(run_cmd[9])
#     beta = float(run_cmd[7])
#     m = float(run_cmd[6])
#     path = run_cmd[10]
#     cwd = os.getcwd()

#     data = pd.read_table(f"{path}/TimeDependence.dat", delimiter="\t", names=["t", "C", "R"])
#     plt.clf()
#     plt.plot(data.R, data.C)
#     plt.xlabel("R")
#     plt.ylabel("C")
#     plt.grid(1)
#     plt.savefig(f"{cwd}/Plots/Phase_alpha{alpha_11:.3f}_beta{beta:.3f}_m{m:.3f}.png")

# def plot_final_states_C1R1(run_cmd):
#     styles = ['r.', 'k.']
#     labels = ['C', None, 'R', None]
#     alpha_11 = float(run_cmd[9])
#     beta = float(run_cmd[7])

#     cwd = os.getcwd()
#     plt.clf()
#     for i in range(len(final_states)):
#         print(i)
#         print(final_states[i][:])
#         plt.plot(mortality_tab, final_states[i][:], styles[int(i/2)], label = labels[i])

#     plt.legend()
#     plt.xlabel("m")
#     plt.ylabel("R,C")
#     plt.grid(1)
#     plt.savefig(f"{cwd}/Plots/Final_state_alpha{alpha_11:.3f}_beta{beta:.3f}.png")



# def main():
#     dt = 5e-2
#     t_max = 1e5
#     C0 = 0.4
#     R0 = 0.7
#     p = 3.
#     #m = .2
#     #beta = 1.
#     b = 2.
#     #alpha_11 = 1.
#     for alpha_11 in alpha_tab:
#         for beta in beta_tab:
#             m_pos = 0
#             for m in mortality_tab:
#                 current_workdir = os.getcwd()
#                 run_folder = f"RUN_alpha{alpha_11:.3f}_beta{beta:.3f}_m{m:.3f}"
#                 output_data = "Output_data"
#                 path = os.path.join(current_workdir, run_folder)
#                 if not os.path.exists(path):
#                     os.mkdir(path) 

#                 run_with_arguments = [current_workdir + '/a.out'] + [str(dt), str(t_max), str(C0),\
#                                                     str(R0), str(p), str(m), str(beta), str(b), str(alpha_11), path]
#                 print(run_with_arguments)
#                 simulation = subprocess.run(run_with_arguments)
#                 plot_phase_C1R1(run_with_arguments)
#                 get_final_states(run_with_arguments, m_pos)
#                 m_pos += 1
#             print(final_states)
#             #plot_final_states_C1R1(run_with_arguments)
