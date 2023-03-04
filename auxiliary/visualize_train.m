function visualize_train(trainset, plot_vars)

    n = size(trainset.X{end},1); %number of variables
    dt = trainset.t{1}(2) - trainset.t{1}(1);
    T = trainset.t{1}(end);
 
    % visualize the predictions for the identfied Koopman model
    load("autokoopman_model.mat", "A","B")

    figure; hold on; box on;
    for r = 1:length(trainset.X)
        %plot Autokoopman vs real trajectory for all simulations
        x = sim_autokoopman(trainset.X{r}(:,1), trainset.XU{r}, @(x) autokoopman(x), A, B, (T/dt)+1);
        if ~any(size(plot_vars)>[1,1]) %singular plot var, plot against time
            plot(trainset.t{r},trainset.X{r}(plot_vars,:),'r:',LineWidth=2);
            plot(trainset.t{r},x(plot_vars,:),'b');
        else
            plot(trainset.X{r}(plot_vars(1),:),trainset.X{r}(plot_vars(2),:),'r:',LineWidth=2);
            plot(x(plot_vars(1),:),x(plot_vars(2),:),'b');
        end
        l = legend('real_trajectory','Autokoopman');
        l.Location = 'best';

    end

end

% Auxiliary Functions -----------------------------------------------------

function x = sim_autokoopman(x0,u,g,A,B,steps)
    x = zeros(size(x0,1),steps);
    x(:,1) = x0;
    for i = 2:steps
        if isempty(B)
            xtemp = A*g(x(:,i-1));
        else
            xtemp = A*g(x(:,i-1)) + B*u(:,i-1);
        end
        x(:,i) = xtemp(1:size(x0,1),:);
    end
end