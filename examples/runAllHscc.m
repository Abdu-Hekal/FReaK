%run all experiments in the paper

%run arch benchmarks
repArch22() %instance 1
repArch22_2() %instance 2

%run aircraft falsification
falsifyAircraft()

% run experiments for different hyper-parameters 
nobservables()
nresets()
timestep()
enhancement()

