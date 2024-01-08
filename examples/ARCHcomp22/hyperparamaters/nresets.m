function nresets()
% nresets - runs experiments in paper on number of resets
%
% Syntax:
%   results = nresets()
%
% Inputs:
%    -
%
% Outputs:
%    results -
%

% Author:       Abdelrahman Hekal
% Written:      30-Oct-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
benches = {}; %empty cell to store benchmarks
%Model transmission benchmark
bench.kfModel = @modelAutoTransmission;
x = stl('x',3);
bench.requirements = {; ...
    "AT1", 1, globally(x(1) < 120,interval(0,20)); ...
    "AT2", 1,globally(x(2) < 4750,interval(0,10)); ...
    "AT6a", 1, implies(globally(x(2)<3000,interval(0,30)),globally(x(1)<35,interval(0,4))); ...
    "AT6b", 1, implies(globally(x(2)<3000,interval(0,30)),globally(x(1)<50,interval(0,8))); ...
    "AT6c", 1, implies(globally(x(2)<3000,interval(0,30)),globally(x(1)<65,interval(0,20))); ...
    "AT6abc", 1, implies(globally(x(2)<3000,interval(0,30)),globally(x(1)<35,interval(0,4))) & implies(globally(x(2)<3000,interval(0,30)),globally(x(1)<50,interval(0,8))) & implies(globally(x(2)<3000,interval(0,30)),globally(x(1)<65,interval(0,20)))
    };
benches{end+1} = bench;
% Chasing cars benchmark
bench.kfModel = @modelCars;
x = stl('x',5);
bench.requirements = {; ...
    "CC1", 10, globally(x(5)-x(4)<=40,interval(0,100)); ...
    "CC2", 10, globally(finally(x(5)-x(4)>=15,interval(0,30)),interval(0,70)); ...
    "CC3", 10, globally(globally(x(2)-x(1)<=20,interval(0,20)) | finally(x(5)-x(4)>=40,interval(0,20)),interval(0,80)); ...
    };
benches{end+1} = bench;

resets=[2,3,5,10,20];
for r = 1:numel(resets)
    solns=dictionary(string.empty,cell.empty);
    for b = 1:length(benches)
        bench = benches{b};
        req = bench.requirements;
        for i = 1:size(req, 1)
            % initialize seeds
            rng(0)
            pyrunfile("seed.py")
            diary('nresets.txt');
            disp("--------------------------------------------------------")
            name = req{i, 1};
            fprintf('Benchmark: %s\n', name);
            fprintf('Number of resets=%d \n',resets(r));
            diary off;

            %initialize progress bar
            msg = sprintf('Runs completed: 0/10');
            fprintf(msg);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));
            for j = 1:10
                kfModel = bench.kfModel();
                kfModel.ak.dt = req{i, 2};
                eq = req{i, 3};
                %settings
                kfModel.nResets=resets(r);

                if name == "NNx"
                    kfModel.U = interval(1.95,2.05);
                end
                if name == "AFC33"
                    kfModel.U = interval([61.2;900],[81.2;1100]);
                end

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
                msg = sprintf('Runs completed: %d/10',j); %Don't forget this semicolon
                fprintf([reverseStr, msg]);
                reverseStr = repmat(sprintf('\b'), 1, length(msg));
            end
            %print info
            fprintf(reverseStr) %remove progress bar
            % Start recording the command line output to a file
            diary('nresets.txt');
            if ~isempty(solns(name))
                printInfo(solns(name),j)
            else
                fprintf('Number of successful falsified traces: 0/%d\n',j)
            end
            % Stop recording the command line output
            diary off;
        end
    end
end

end

