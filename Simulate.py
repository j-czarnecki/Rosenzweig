import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os 
import subprocess

N_EQUATIONS = 2
alpha_tab = [0.5, 0.6, 0.7, 0.8, 0.9, 1.]
beta_tab = [0.2, 0.4, 0.6, 0.8, 0.9, 1.]
#beta_tab = [1.]
mortality_tab = [0.05*i for i in range(1,21)]
final_states = [[0]*len(mortality_tab) for i in range(2*N_EQUATIONS)]

# alpha_tab = [1.]
# beta_tab = [1.]
fig = plt.figure()



def plot_phase_C1R1(run_cmd):
    alpha_11 = float(run_cmd[9])
    beta = float(run_cmd[7])
    m = float(run_cmd[6])
    path = run_cmd[10]
    cwd = os.getcwd()

    data = pd.read_table(f"{path}/TimeDependence.dat", delimiter="\t", names=["t", "C", "R"])
    plt.clf()
    plt.plot(data.R, data.C)
    plt.xlabel("R")
    plt.ylabel("C")
    plt.grid(1)
    plt.savefig(f"{cwd}/Plots/Phase_alpha{alpha_11:.3f}_beta{beta:.3f}_m{m:.3f}.png")

def plot_final_states_C1R1(run_cmd):
    styles = ['r.', 'k.']
    labels = ['C', None, 'R', None]
    alpha_11 = float(run_cmd[9])
    beta = float(run_cmd[7])

    cwd = os.getcwd()
    plt.clf()
    for i in range(len(final_states)):
        print(i)
        print(final_states[i][:])
        plt.plot(mortality_tab, final_states[i][:], styles[int(i/2)], label = labels[i])

    plt.legend()
    plt.xlabel("m")
    plt.ylabel("R,C")
    plt.grid(1)
    plt.savefig(f"{cwd}/Plots/Final_state_alpha{alpha_11:.3f}_beta{beta:.3f}.png")



def get_final_states(run_cmd, m_pos):
    path = run_cmd[10]
    alpha_11 = float(run_cmd[9])
    beta = float(run_cmd[7])
    m = float(run_cmd[6])
    print("Alpha ", alpha_11)
    print("Beta ", beta)
    print("m ", m)

    data = pd.read_table(f"{path}/FinalState.dat", delimiter="\t", names=["maximum", "minimum"])
    for i in range(len(data.maximum)):
        final_states[2*i][m_pos] = data.maximum[i]
        final_states[2*i+1][m_pos] = data.minimum[i]


def main():

    dt = 5e-2
    t_max = 1e5
    C0 = 0.4
    R0 = 0.7
    p = 3.
    #m = 1.
    #beta = 1.
    b = 2.
    #alpha_11 = 1.
    for alpha_11 in alpha_tab:
        for beta in beta_tab:
            m_pos = 0
            for m in mortality_tab:
                current_workdir = os.getcwd()
                run_folder = f"RUN_alpha{alpha_11:.3f}_beta{beta:.3f}_m{m:.3f}"
                output_data = "Output_data"
                path = os.path.join(current_workdir, run_folder)
                if not os.path.exists(path):
                    os.mkdir(path) 

                run_with_arguments = [current_workdir + '/a.out'] + [str(dt), str(t_max), str(C0),\
                                                    str(R0), str(p), str(m), str(beta), str(b), str(alpha_11), path]
                print(run_with_arguments)
                simulation = subprocess.run(run_with_arguments)
                plot_phase_C1R1(run_with_arguments)
                get_final_states(run_with_arguments, m_pos)
                m_pos += 1
            print(final_states)
            #plot_final_states_C1R1(run_with_arguments)

if __name__ == "__main__":
    main()
