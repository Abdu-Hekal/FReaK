function figure_settings(fig)

ax=gca(fig);

set(ax, 'fontname', 'Times New Roman', 'FontSize', 5);

set(fig, 'units', 'centimeters', 'position', [1 1 8.89 6.5]);
set(fig, 'Color', 'w');

ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
ax.ZMinorTick = 'on';

end