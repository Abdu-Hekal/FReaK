function figure_settings(fig)

ax=gca(fig);

set(fig, 'Color', 'w');

%figure settings for HSCC paper 
set(fig, 'units', 'centimeters', 'position', [1 1 8.89 6.5]);
set(ax, 'fontname', 'Times New Roman', 'FontSize', 7);

ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
ax.ZMinorTick = 'on';

%change legend size
L = legend;
L.ItemTokenSize(1) = 10;

end