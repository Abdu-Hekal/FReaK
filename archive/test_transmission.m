dt = 0.01;
figure; hold on; box on;
model = model_AutoTransmission();
x0 = (model.R0.sup-model.R0.inf)*rand()+model.R0.inf;
u = [linspace(0,model.T-model.T/model.N,model.N)',(zeros(1,model.N))'+100,(zeros(1,model.N))'+325]; %check this
[t, x] = run_simulation(model.name, model.T, x0, u);
plot(x(:,1),x(:,2),'r');

% visualize the predictions for the identfied Koopman model
load("autokoopman_model.mat", "A","B")
x = sim_autokoopman(x0,u',@(x) autokoopman(x), A, B, model.T/dt);
plot(x(1,:),x(2,:),'b');

l = legend('real_trajectory','Autokoopman');
l.Location = 'best';

figure; hold on; box on;
model = model_AutoTransmission();
x0 = (model.R0.sup-model.R0.inf)*rand()+model.R0.inf;
u = [linspace(0,model.T-model.T/model.N,model.N)',(zeros(1,model.N))'+100,(zeros(1,model.N))']; %check this
[t, x] = run_simulation(model.name, model.T, x0, u);
plot(x(:,1),x(:,2),'r');

% visualize the predictions for the identfied Koopman model
load("autokoopman_model.mat", "A","B")
x = sim_autokoopman(x0,u',@(x) autokoopman(x), A, B, model.T/dt);
plot(x(1,:),x(2,:),'b');

l = legend('real_trajectory','Autokoopman');
l.Location = 'best';




function x = sim_autokoopman(x0,u,g,A,B,steps)
    x = zeros(size(x0,1),steps);
    x(:,1) = x0;
    for i = 2:steps
        if isempty(B)
            xtemp = A*g(x(:,i-1));
        else
            xtemp = A*g(x(:,i-1)) + B*u(2:end,i-1);
        end
        x(:,i) = xtemp(1:size(x0,1),:);
    end
end