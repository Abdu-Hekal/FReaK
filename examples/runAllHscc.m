%run all experiments in the paper

%run arch benchmarks
disp('******************* Instance 1 *******************')
repArch22() %instance 1
disp('******************* Instance 2 *******************')
repArch22_2() %instance 2

%run aircraft falsification
disp('******************* Aircraft *******************')
falsifyAircraft()
disp('******************* Aircraft - Staliro *******************')
ah_aircraftODE()

% run experiments for different hyper-parameters 
disp('******************* Number of Observables *******************')
nobservables()
disp('******************* Resets *******************')
nresets()
disp('******************* Timestep *******************')
timestep()
disp('******************* Enhancements *******************')
enhancement()

