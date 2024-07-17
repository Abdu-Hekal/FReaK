bench.kfModel = @modelAircraftODE;
x = stl('x',3);
bench.requirements = {; ...
    "phi1", implies(globally(x(1) >=250 & x(1) <=260,interval(1,1.5)),globally(~(x(1)>=230 & x(1)<=240),interval(3,4))); ...
    "phi2",globally(x(3)>0,interval(0,4)); ...
    };

% Start recording the command line output to a file
diary('aircraft.txt');

solns=dictionary(string.empty,cell.empty);
req = bench.requirements;
for i = 1:size(req, 1)
    % initialize seeds
    rng(0)
    pyrunfile("seed.py")
    disp("--------------------------------------------------------")
    kfModel = bench.kfModel();
    name = req{i, 1};
    eq = req{i, 2};

    kfModel.spec = specification(eq,'logic');
    kfModel.runs=10;
    kfModel.verb=1;
    kfSolns = falsify(kfModel);

    %print info
    fprintf('Benchmark: %s\n', name);
    printInfo(kfSolns)

end
diary off;


