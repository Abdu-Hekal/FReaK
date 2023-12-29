%setup koopman falsification
dockerSetupKF;
%run all experiments
warning('off');
runAllHscc;