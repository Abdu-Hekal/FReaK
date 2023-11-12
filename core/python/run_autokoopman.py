from scipy.io import savemat
from autokoopman import auto_koopman

import autokoopman.core.trajectory as traj
import numpy as np
import torch
import random

def run(times, trajectories, param_dict, inputs_list):

    training_data = []
    for i, (time, trajectory) in enumerate(zip(times, trajectories)):
        inputs = np.atleast_2d(inputs_list[i]).T if inputs_list else None
        training_data.append(traj.Trajectory(np.asarray(time), np.atleast_2d(trajectory).T, inputs))
        
    ids = np.arange(0, len(training_data)).tolist()
    training_data = traj.TrajectoriesData(dict(zip(ids, training_data)))

    if training_data.n_trajs > 3:
        n_splits = int(training_data.n_trajs/2) if training_data.n_trajs%2==0 else None
    else:
        n_splits=None

    experiment_results = auto_koopman(
        training_data,  # list of trajectories
        sampling_period=param_dict["dt"],  # sampling period of trajectory snapshots
        obs_type=param_dict["obsType"],  # use Random Fourier Features Observables
        opt=param_dict["opt"],  # grid search to find best hyperparameters
        n_obs=int(param_dict["nObs"]),  # maximum number of observables to try
        max_opt_iter=200,  # maximum number of optimization iterations
        grid_param_slices=int(param_dict["gridSlices"]),
        # for grid search, number of slices for each parameter
        rank=tuple(param_dict["rank"]), # rank range (start, stop, step) DMD hyperparameter
        n_splits=n_splits,
        verbose= False,
    )
    
    
    model = experiment_results['tuned_model']
    # get evolution matrices
    A, B = model.A, model.B
    w = model.obs_func.observables[1].w
    u = model.obs_func.observables[1].u
    
    koopman_model = {"A": A, "B": B, "w": w, "u": u}
    savemat("autokoopman_model.mat", koopman_model)

    params = experiment_results['hyperparameters']
    paramVals = experiment_results['hyperparameter_values']

#     print(params)
#     print(paramVals)

    return koopman_model


koopman_model = run(times,trajectories, param_dict,inputs_list)



