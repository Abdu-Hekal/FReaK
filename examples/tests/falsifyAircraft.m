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
    %name and stl
    name = req{i, 1};
    eq = req{i, 2};
    disp("--------------------------------------------------------")
    fprintf('Benchmark: %s\n', name);
    %initialize progress bar
    msg = sprintf('Runs completed: 0/10');
    fprintf(msg);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    for j = 1:10
        kfModel = bench.kfModel();

        kfModel.spec = specification(eq,'logic');
        kfSoln = falsify(kfModel);

        if j==1
            solns(name)={{}};
        end
        if kfSoln.falsified
            soln=solns(name);
            soln{1}{end+1}=kfSoln;
            solns(name)=soln;
        end
        % Display the progress
        msg = sprintf('Runs completed: %d/10',j);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end
    %print info
    fprintf(reverseStr) %remove progress bar
    if ~isempty(solns(name))
        printInfo(solns(name),j)
    else
        fprintf('Number of successful falsified traces: 0/%d\n',j)
    end
end
diary off;


