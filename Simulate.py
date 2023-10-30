import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import pandas as pd
import numpy as np
import os 
import subprocess
import matplotlib.collections as mcoll
import matplotlib.path as mpath
import math


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

def main():
    dt = 1e-1
    t_max = 2e6
    C0 = 0.4
    R1_0 = 0.71
    R2_0 = 0.7 
    beta0 = 0.49
    gamma0 = 0.
    p = 3.
    m = 0.2
    b = 2.
    alpha_11 = 0.95
    alpha_12 = 0.065
    alpha_21 = 0.065
    alpha_22 = 0.95
    v = 80e-8

    for alpha_11 in alpha_tab:
        alpha_22 = alpha_11
        v_pos = 0
        for v in v_tab:
            current_workdir = os.getcwd()
            path = f"/mnt/HDD1TB/Rosenzweig/RUN_alpha{alpha_11:.3f}_v{v*1e8:.3f}"
            if not os.path.exists(path):
                os.mkdir(path) 

            run_with_arguments = [current_workdir + '/a.out'] + [str(dt), str(t_max), str(C0),\
                                            str(R1_0), str(R2_0), str(beta0), str(gamma0), str(p), str(m), str(b), \
                                            str(alpha_11), str(alpha_12), str(alpha_21), str(alpha_22), str(v), path]
            print(run_with_arguments)
            simulation = subprocess.run(run_with_arguments)
            plot_time_dependence(run_with_arguments)
            plot_phase(run_with_arguments)
            #get_final_states(run_with_arguments, v_pos)
            get_final_phase(run_with_arguments, v_pos)
            v_pos += 1
        #plot_final_states(run_with_arguments)
        #clear_final_states()
        plot_final_phase(run_with_arguments)
        clear_final_phase()


if __name__ == "__main__":
    main()




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
