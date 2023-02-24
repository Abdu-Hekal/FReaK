from scipy.io import savemat
from autokoopman import auto_koopman

import autokoopman.core.trajectory as traj
import numpy as np

def run(times, trajectories, param_dict, inputs_list):

    data = []
    for i, (time, trajectory) in enumerate(zip(times, trajectories)):
        inputs = np.asarray(inputs_list[0]).T if inputs_list else None
        data.append(traj.Trajectory(np.asarray(time), np.asarray(trajectory).T, inputs))
    
    if param_dict["obs_type"] == 'deep':
        opt = 'bopt'
    else:
        opt = 'grid'
    
    
    training_data = [data[i] for i in range(0, len(data))]
    
    ids = np.arange(0, len(training_data)).tolist()
    training_data = traj.TrajectoriesData(dict(zip(ids, training_data)))

    experiment_results = auto_koopman(
        training_data,  # list of trajectories
        sampling_period=param_dict["samp_period"],  # sampling period of trajectory snapshots
        obs_type=param_dict["obs_type"],  # use Random Fourier Features Observables
        opt=opt,  # grid search to find best hyperparameters
        n_obs=int(param_dict["n_obs"]),  # maximum number of observables to try
        max_opt_iter=200,  # maximum number of optimization iterations
        grid_param_slices=int(param_dict["grid_param_slices"]),
        # for grid search, number of slices for each parameter
        rank=(1, 20, 1) # rank range (start, stop, step) DMD hyperparameter
    )
    
    
    model = experiment_results['tuned_model']
    # get evolution matrices
    A, B = model.A, model.B
    w = model.obs_func.observables[1].w
    u = model.obs_func.observables[1].u
    
    koopman_model = {"A": A, "B": B, "w": w, "u": u}
    savemat("autokoopman_model.mat", koopman_model)

    trajectory = model.solve_ivp(
    initial_state=training_data.get_trajectory(0).states[0],
    tspan=(0.0, 7.0),
    sampling_period=0.01
    )
    print(trajectory.states)

    return koopman_model

koopman_model = run(times,trajectories, param_dict,inputs_list)



