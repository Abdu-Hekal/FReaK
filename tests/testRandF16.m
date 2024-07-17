kf=modelF16;
kf.verb=2;
kf.maxSims=100000;

[soln,sims]=kf.randFalsify;
save('soln.mat','soln')
save('sims.mat','sims')

X0=[];
for ii=1:numel(sims.X)
    x=sims.X{ii};
    X0=[X0;x(1,:)];
end
x=X0(:,4); y=X0(:,5); z=X0(:,6);
scSize = max(min(ceil(2000/size(x,1)),100),10);
%create 3d scatter plot
figure;
scatter3(x, y, z, scSize, cell2mat(sims.ROB)', 'filled');
% Customize plot
xlabel('Var 4');
ylabel('Var 5');
zlabel('Var 6');
title('Heatmap for Rand Falsification of F16');
colormap(flip(jet));  %colormap
colorbar;

% Create 2s scatter plots with color encoding
% Create a figure
figure;

% Scatter plot for (x, y)
subplot(2, 2, 1);
scatter(x, y, scSize, cell2mat(sims.ROB)', 'filled');
title('(x, y)');
xlabel('X-Axis');
ylabel('Y-Axis');
% Apply the same colormap
colormap(flip(jet));
colorbar('Location', 'eastoutside');

% Scatter plot for (x, z)
subplot(2, 2, 2);
scatter(x, z, scSize, cell2mat(sims.ROB)', 'filled');
title('(x, z)');
xlabel('X-Axis');
ylabel('Z-Axis');
% Apply the same colormap
colormap(flip(jet));
colorbar('Location', 'eastoutside');

% Scatter plot for (y, z)
subplot(2, 2, 3);
scatter(y, z, scSize, cell2mat(sims.ROB)', 'filled');
title('(y, z)');
xlabel('Y-Axis');
ylabel('Z-Axis');
% Apply the same colormap
colormap(flip(jet));
colorbar('Location', 'eastoutside');


% Adjust layout
sgtitle('All Combinations of (x, y), (x, z), and (y, z) in 2D Scatter Plots');

% Find the indices of the minimum value
[~, minIndex] = min(cell2mat(sims.ROB));
minX = x(minIndex)
minY = y(minIndex)
minZ = z(minIndex)

