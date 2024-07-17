from scipy.io import savemat
from autokoopman import auto_koopman
from autokoopman.core.format import hide_prints
from sklearn.preprocessing import normalize, MinMaxScaler, RobustScaler

import autokoopman.core.trajectory as traj
import numpy as np
import random

def run(times, trajectories, param_dict, inputs_list, rob_list,state_weights):

    # Sort all lists based on robustness
#     if inputs_list:
#         sorted_lists = sorted(zip(rob_list, trajectories, times,inputs_list), key=lambda x: x[0])
#         rob_list, trajectories, times, inputs_list = zip(*sorted_lists)
#     else:
#         sorted_lists = sorted(zip(rob_list, trajectories, times), key=lambda x: x[0])
#         rob_list, trajectories, times = zip(*sorted_lists)

    if rob_list is not None:
        rob_list=np.array(rob_list).reshape(-1,1)
        state_weights=np.array(state_weights)
        inputs_list=np.array(inputs_list)
        cost_func='weighted'
#         # remove large outliers (large +ve robustness)
#         q1 = np.percentile(rob_list, 25)
#         q3 = np.percentile(rob_list, 75)
#         iqr = q3 - q1
#         upper_bound = q3 + 1.5 * iqr
#         indices=np.where(rob_list<=upper_bound)[0]
#         rob_list=rob_list[indices]
#         trajectories=np.array(trajectories)[indices]
#         times=np.array(times)[indices]
#         inputs_list=np.array(inputs_list)
#         if inputs_list.any():
#             inputs_list=inputs_list[indices]
# #         minmax normalize robustness
#         scaler = MinMaxScaler()
#         rob_list = scaler.fit_transform(rob_list.reshape(-1,1))
    else:
        cost_func='total'
        inputs_list=np.array(inputs_list)

    training_data = []
    weights = []
    for i, (time, trajectory) in enumerate(zip(times, trajectories)):
        # note we are currently not using time list as we learn a uniform time trajectory, however we may want to change this if we learn a continous model or other.
        inputs = np.atleast_2d(inputs_list[i]).T if inputs_list.any() else None
        training_traj=traj.UniformTimeTrajectory(np.atleast_2d(trajectory).T, inputs, param_dict["dt"])
        training_data.append(training_traj)
        if rob_list is not None:
            w = np.ones(training_traj.states.shape)*state_weights*(1/rob_list[i]) #(i+1/len(trajectories))/(rob_list[i]+(len(trajectories)))
            weights.append(w)

    #convert to np array for data manipulation
    training_data=np.array(training_data)
    weights=np.array(weights)
    
    if rob_list is None:
        weights=None

    ids = np.arange(0, len(training_data)).tolist()
    training_data = traj.UniformTimeTrajectoriesData(dict(zip(ids, training_data)))
    if rob_list is not None:
        weights=dict(zip(ids, weights))

    if training_data.n_trajs > 3:
#         n_splits = int(training_data.n_trajs/3) if training_data.n_trajs%3==0 else None
        n_splits = int(training_data.n_trajs/2) if training_data.n_trajs%2==0 else None
    else:
        n_splits=None
   
    experiment_results = auto_koopman(
        training_data,  # list of trajectories
        sampling_period=param_dict["dt"],  # sampling period of trajectory snapshots
        obs_type=param_dict["obsType"],  # use Random Fourier Features Observables
        cost_func=cost_func,   # use "weighted" cost function
        learning_weights=weights, # weight the eDMD algorithm objectives
        scoring_weights=weights, # pass weights as required for cost_func="weighted"
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


koopman_model = run(times,trajectories, param_dict,inputs_list,rob_list,state_weights)