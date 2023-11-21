load("autokoopman_model.mat", "A","B")
g = @(x) autokoopman(x);

model = @(t,x,u) [A*x(1:size(A,2))+B*u(1:size(B,2))];
disp(model(1,g([0;1000;1]),[1;2]))
options = staliro_options();
options.spec_space='X';

phi = '!pred1';
preds(1).str='pred1';
preds(1).A = [-1,zeros(1,size(A,2)-1)];
size(preds(1).A)
preds(1).b=[-1.5];
run = staliro(model, [g([0;1000;1]),g([0;1000;1])], [0 100;0 325], [3000;3000], phi, preds, 10, options);
run.bestSample

function [ret] = create_lin_model(t,x,u)
load("autokoopman_model.mat", "A","B")
g = @(x) autokoopman(x);

f=[];
for row=1:length(A)
    f = [f,A(row)*x+B(row)*u];
end
end