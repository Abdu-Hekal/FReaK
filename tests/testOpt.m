% Define signal temporal logic formula
x = stl('x',4);
stl_formula = globally(x(1)<0.5 & (x(2)>0.3 | x(3)<0),interval(0,1));

% Define time vector
t = 0:0.1:1;

% Define linear constraints
A = [1 0; 0 1; -1 0; 0 -1];
b = [1; 1; 0; 0];
constraints = struct('A',A,'b',b);

% Solve MILP
[status, x, obj_val] = milp_robustness(stl_formula, t, constraints);


% Print results
if status == 1
    for i=1:length(x)
        fprintf('Optimal input signal: x%d =%f\n', i, x(i));
    end
end
