function corePlotReach(critU,critX,t,xt,x0,A,B,g,R)
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
    p1=plot(x(plotVars(1),1:end),x(plotVars(2),1:end),'Color',PS.DGrey5,'LineWidth',1,'DisplayName','Koopman trajectory');
    p2=plot(critX(:,plotVars(1)),critX(:,plotVars(2)), 'Color' , PS.MyGreen3,'LineWidth',1,'DisplayName','Real trajectory');

    for i=1:size(xt)-1
        if i==1
            plot(R.zono{i},plotVars,'FaceColor',PS.MyBlue3,'FaceAlpha',.02,'DisplayName','Reachable set')
        else
            plot(R.zono{i},plotVars,'FaceColor',PS.MyBlue3,'FaceAlpha',.02,'HandleVisibility','off')
        end
    end

    %custom limits
    xlim([-10, 80])
    ylim([-1000, 5000])

    % Find the index of the maximum value in y
    [~, maxIndex] = max(critX(1:end,plotVars(2)));
    % Get the corresponding x and y values for the maximum point
    xMax = critX(maxIndex,plotVars(1));
    yMax = critX(maxIndex,plotVars(2));
    % Draw a horizontal line at ymax
    line([min(xlim), xMax], [yMax, yMax],'Color',PS.MyBlack,'LineStyle', '--','HandleVisibility','off');

    % Add a custom y-axis value at the maximum y without modifying existing y-ticks
    yticks = get(gca, 'YTick');
    yticks = [yticks,round(yMax)];
    yticks=sort(yticks);
    set(gca, 'YTick', yticks);

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


    export_fig reach.pdf


    drawnow;
end
end