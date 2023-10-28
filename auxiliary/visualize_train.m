function visualize_train(varargin)
%varargin: trainset,plot_vars,xlabel,ylabel

numArgs = length(varargin);
assert(numArgs>=2,'must have at least 2 args, trainset and plot_vars')
trainset=varargin{1};
plot_vars=varargin{2};

n = size(trainset.X{end},1); %number of variables
dt = trainset.t{1}(2) - trainset.t{1}(1);
T = trainset.t{1}(end);

% visualize the predictions for the identfied Koopman model
load("autokoopman_model.mat", "A","B")

fig=figure; hold on; box on;

% The standard values for colors saved in PLOT_STANDARDS() will be accessed from the variable PS
PS = PLOT_STANDARDS();
%settings for figure
figure_settings(fig);

for r = 1:length(trainset.X)
    %plot Autokoopman vs real trajectory for all simulations
    x = sim_autokoopman(trainset.X{r}(:,1), trainset.XU{r}, @(x) autokoopman(x), A, B, (T/dt)+1);

    if ~any(size(plot_vars)>[1,1]) %singular plot var, plot against time
        p1=plot(trainset.t{r},x(plot_vars,:));
        p2=plot(trainset.t{r},trainset.X{r}(plot_vars,:));

    else
        p1=plot(x(plot_vars(1),:),x(plot_vars(2),:));

        p2=plot(trainset.X{r}(plot_vars(1),:),trainset.X{r}(plot_vars(2),:));
    end
    %style plots 
    set(p1, 'LineWidth', 1, 'Color', PS.DGrey5);
    set(p2, 'LineWidth', 1, 'Color',PS.MyGreen3);


    %add legend
    l = legend('Koopman trajectory','Real trajectory');
    l.Location = 'northoutside';
    l.NumColumns = 2;

end
% Add axis titles
if numArgs>2
    xlabel(varargin{3});
end
if numArgs>3
    ylabel(varargin{4});
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