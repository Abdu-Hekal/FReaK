%Ensure that autokoopman is installed & imported in your python environment
py.importlib.import_module('autokoopman');

% initalize training data
X = {}; U={}; T=7; t = {}; max_train_size=20; 

%generate random initial set
x = 0.3*rand()+1.25;
y = 0.1*rand()+2.25;
init_set = [x;y];

for i = 1:max_train_size
    [sys, X, t, x, crit_x] = vanderpolSymbolicRFF(init_set,X, t, T);
    %retrain with initial set as the critical set found in prev iteration
    init_set=crit_x(1,:)';
    disp(i)
    if any(crit_x(:,1)>=2.095)
        break;
    end
end

train_iter = ['training iterations required: ',num2str(i)];
disp(train_iter)

%plot falsifying trace
figure; hold on; box on;
plot([2.095 2.095],[-3, 3],'--r')
plot(crit_x(:,1),crit_x(:,2),'g');

figure; hold on; box on;
for r = 1:length(X)
    %plot Autokoopman vs real trajectory for all simulations
    U = zeros(size(X{r})-1);
    x = simulateKoopman(sys,@(x) autokoopmanVanderPol(x),X{r}(:,1),U);
    plot(x(1,:),x(2,:),'b');
    plot(X{r}(1,:),X{r}(2,:),'r');

    legend('Autokoopman','real_trajectory')

end