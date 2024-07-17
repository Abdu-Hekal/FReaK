function corePlotFalsify(critU,critX,t,xt,x0,A,B,g,R)
fig=figure;
% The standard values for colors saved in PLOT_STANDARDS() will be accessed from the variable PS
PS = PLOT_STANDARDS();
%settings for figure
figure_settings(fig);
set(fig, 'units', 'centimeters', 'position', [1 1 8.89 8.5]);


plotVars=[1,2]; %[3];
drawu=critU(:,2:end)';
x = g(x0);
for i = 1:size(xt)-1
    if ~isempty(drawu)
        x = [x, A*x(:,end) + B*drawu(:,i)];
    else
        x = [x, A*x(:,end)];
    end
end
hold on; box on;
if ~any(size(plotVars)>[1,1]) %singular plot var, plot against time
    plot(xt,x(plotVars(1),1:end),'r','LineWidth',2);
    plot(t,critX(1:end,plotVars(1)),'g','LineWidth',2)
else
    p1=plot(x(plotVars(1),1:end),x(plotVars(2),1:end),'Color',PS.DGrey5,'LineWidth',1,'DisplayName','Koopman traj');
    from=1; to=400;
    p2=plot(critX(from:to,plotVars(1)),critX(from:to,plotVars(2)),'Color' ,PS.MyGreen3,'LineWidth',1,'DisplayName','Falsifying traj (t\leq4)');
    p3=plot(critX(to:end,plotVars(1)),critX(to:end,plotVars(2)), 'Color' , PS.MyGreen3,'LineStyle','--','LineWidth',1,'DisplayName','Falsifying traj (t\geq4)');

    %custom limits
    xlim([-10, 80])
    ylim([-1000, 5000])


    % plot unqualified section
    % Get the current x-limits of the plot
    xLimits = xlim;
    % Define the boundary of the half-space (y > 3000)
    y_boundary = 3000;
    % Plot the half-space based on the current x-limits
    x = linspace(xLimits(1), xLimits(2), 100); % Automatically selected x-range
    y = y_boundary * ones(size(x)); % y = 3000 for all x in the selected range
    %plot halfspace
    fill([x, fliplr(x)], [y, max(ylim)*ones(size(y))], PS.DGrey5, 'FaceAlpha', 0.15,'DisplayName','Invalid region');

    xlabel('Speed');
    ylabel('Angular velocity');

    %add legend
    l = legend;
    l.Location = 'northoutside';
    l.NumColumns = 3;

    % plot unsafe section
    % Get the current y-limits of the plot
    yLimits = ylim;
    yLimits(2) = 3000;
    % Define the boundary of the half-space (x > 35)
    x_boundary = 35;
    % Plot the half-space based on the current x-limits
    y = linspace(yLimits(1), yLimits(2), 100); % Automatically selected y-range
    x = x_boundary * ones(size(x)); % x = 35 for all x in the selected range
    %plot halfspace
    fill([x, max(xlim)*ones(size(x))], [y, fliplr(y)], 'r', 'FaceAlpha', 0.05,'DisplayName','Unsafe region');

    %bring trajectories to front
    uistack(p1, 'top');
    uistack(p2, 'top');
    uistack(p3, 'top');

    export_fig falsification.pdf


    drawnow;
end
end