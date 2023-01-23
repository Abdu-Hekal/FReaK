% load your python environment that autokoopman is installed & imported in
pyenv("Version",'/Users/b6062805/Documents/Koopman/autokoopman_vitualenv/bin/python','ExecutionMode','InProcess');
py.importlib.import_module('autokoopman');

% initalize training data
X = {}; U={}; T=7; t = {};

[sys, X, U, t, x, crit_x] = vanderpolSymbolicRFF(X, U, t, T);
