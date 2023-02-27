figure; hold on; box on;
model = model_AutoTransmission();
for i = 1:1
    x0 = (model.R0.sup-model.R0.inf)*rand()+model.R0.inf;
    u = [linspace(0,model.T-model.T/model.N,model.N)',(zeros(1,model.N))'+100,(zeros(1,model.N))']; %check this
    [t, x] = run_simulation(model.name, model.T, x0, u);
    plot(x(:,1),x(:,2),'b');
end