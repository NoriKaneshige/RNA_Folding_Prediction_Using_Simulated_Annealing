import numpy as np
import matplotlib.pyplot as plt
# ===============================================
# TODO: you can import more functions from PyRNA as needed
from PyRNA import perturb_stem, get_free_energy_from_stems, stem2basepair_matrix, visualize_structure
# ===============================================
import copy


def muller_potential(x, y, openMM=True):
    """ Defines the 2D Muller potential """
    if openMM:
        # https://gist.github.com/rmcgibbo/6094172
        a = np.array([-1, -1, -6.5, 0.7])
        b = np.array([0, 0, 11, 0.6])
        c = np.array([-10, -10, -6.5, 0.7])
        A = np.array([-200, -100, -170, 15])
        x0 = np.array([1, 0, -0.5, -1])
        y0 = np.array([0, 0.5, 1.5, 1])
    else:
        A = np.array([0, -10, -10, 0.5])
        a = np.array([-1, -1, -6.5, 0.7])
        b = np.array([0, 0, 11, 0.6])
        c = np.array([-10, -10, -6.5, 0.7])
        x0 = np.array([1, 0, -0.5, -1])
        y0 = np.array([0, 0.5, 1.5, 1])

    # accumulate potentials
    return(np.sum(A*np.exp(a*(x-x0)**2+b*(x-x0)*(y-y0)+c*(y-y0)**2)))


def muller_potential_grid(X, Y, openMM=True):
    """ Computes the 2D Muller potential over a grid """
    # get shape of the grid and then compute and store results
    nr, nc = X.shape[0], X.shape[1]
    potential = np.zeros((nr, nc))
    for i in range(nr):
        for j in range(nc):
            x, y = X[i, j], Y[i, j]
            potential[i, j] = muller_potential(x, y, openMM=True)

    # return the potential grid
    return(potential)


def plot_muller_potential(x_min=-2.0, x_max=1.0, y_min=-1.0, y_max=2.5, delta=0.025, origin='lower', colors=plt.cm.viridis):
    """ Plots the 2D Muller potential """
    # set up grid and compute potential
    x = np.arange(x_min, x_max, delta)
    y = np.arange(y_min, y_max, delta)
    X, Y = np.meshgrid(x, y)
    Z = muller_potential_grid(X, Y)
    nr, nc = Z.shape

    # plot contours
    fig1, ax2 = plt.subplots(constrained_layout=True)
    CS = ax2.contourf(X, Y, Z.clip(max=120), 10, cmap=colors, origin=origin)
    CS2 = ax2.contour(CS, levels=CS.levels[::1], colors='r', origin=origin)

    # labels
    ax2.set_title('Muller Potential')
    ax2.set_xlabel('Coordinate-X')
    ax2.set_ylabel('Coordinate-Y')

    # Make a colorbar for the ContourSet returned by the contourf call.
    cbar = fig1.colorbar(CS)
    cbar.ax.set_ylabel('Energy')

    # Add the contour line levels to the colorbar
    cbar.add_lines(CS2)


def energy(state, system="2D_muller"):
    """ Returns the energy of the system """
    if system == "2D_muller":
        energy = muller_potential(state[0], state[1])
    if system == "RNA":
        # ===============================================
        # TODO: replace the pass with your code

        energy = get_free_energy_from_stems(state['assembled_stems'], state['stem_energies'])
        
        """ state['stems'] just returns the list of all stems. This next changes, 
        so every time you call this function, with these arguments, you get the same energies!
        You need to find the list that stores the stems in the current state of the RNA."""

        # ===============================================
    return(energy)


def perturb_state(state, system="2D_muller", step_size=0.5):
    """ Randomly perturbs the state of a system """
    if system == "2D_muller":
        state_tmp = perturb_position(state, step_size)
    if system == "RNA":
        # ===============================================
        # TODO: replace the pass with your code

        state_tmp = perturb_stem(state)

        # ===============================================
    return(state_tmp)


def perturb_position(state, step_size=0.5):
    " Randomly perturbs the position of a 2D object "
    state[0] += step_size*np.random.uniform(-1, 1)
    state[1] += step_size*np.random.uniform(-1, 1)
    return(state)


def Metropolis_criteria(state0, state1, temperature, distribution_parameter, system='2D_muller'):
    """ Executes the Metropolis_criteria, i.e., decides whether to accept or reject a move """
    r = np.random.uniform(0, 1)
    energy0 = energy(state0, system)
    energy1 = energy(state1, system)
    delta_energy = energy1 - energy0
    prob = np.exp(-delta_energy/(temperature * distribution_parameter))
    if r < prob or delta_energy < 0.0:
        state = state1
    else:
        state = state0
    return(state)


def update_distribution_parameter(distribution_parameter, cooling_rate):
    """ Updates the simulated annealing distribution parameter """
    return(distribution_parameter*cooling_rate)


def print_output(i, j, state0, state1, distribution_parameter0, system='2D_muller'):
    """ Print the state information out to screen """
    if system == "2D_muller":
        print("Position at iteration=%s and Cooling step=%s: X=%4.2f Y=%4.2f Energy=%4.2f/%4.2f Distribution Parameter=%4.5f" %
              (i, j, state0[0], state0[1], energy(state0, system), energy(state1, system), distribution_parameter0))
    if system == "RNA":
        # visualize structure
        label = "Fold at iteration=%s and Cooling step=%s: Energy=%4.2f/%4.2f Distribution Parameter=%4.5f" % (
            i, j, energy(state0, system), energy(state1, system), distribution_parameter0)
        matrix = stem2basepair_matrix(
            state0['sequence'], state0['assembled_stems'], state0['stems_s1'], state0['stems_s2'])
        visualize_structure(matrix, label)


def simulated_annealing(initial_state=[0, 0],
                        temperature=1,
                        distribution_parameter=100,
                        step_size=0.5,
                        iterations=100,
                        equilibration_steps=1000,
                        cooling_steps=100,
                        cooling_rate=0.95, debug=False, system="2D_muller"):
    """ Runs simulated annealing """
    # initialize state
    state0 = copy.deepcopy(initial_state)
    
    current_state_list = []
    current_state_energy_list = []
    # Monte Carlo Loop
    # ===============================================
    # TODO: Modify this function following instructions in "RNA_folding_assignment.ipynb"
    # ===============================================
    for i in range(iterations):

        current_state = None
        # initialize distribution_parameter
        distribution_parameter0 = copy.deepcopy(distribution_parameter)
        for j in range(cooling_steps):
            for k in range(equilibration_steps):
                # perturb state
                state1 = perturb_state(copy.deepcopy(state0), system, step_size)

                # update state (if accepted)
                state0 = Metropolis_criteria(copy.deepcopy(state0), copy.deepcopy(
                    state1), temperature, distribution_parameter0, system)
                current_state = state0
                # print(current_state["assembled_stems"])
                if debug:
                    print_output(i, j, state0, state1, distribution_parameter0, system)

                # update the distribution parameter
                distribution_parameter0 = update_distribution_parameter(
                    distribution_parameter0, cooling_rate)

        # Get final position and energy
        if not debug:
            print_output(i, j, state0, state1, distribution_parameter0, system)

        current_state_list.append(current_state)
        current_state_energy_list.append(energy(current_state, system="RNA"))
        print("assembled_stems: {}, current_state_energies: {}".format(current_state['assembled_stems'], energy(current_state, system="RNA")))

    return(current_state_list,current_state_energy_list)
