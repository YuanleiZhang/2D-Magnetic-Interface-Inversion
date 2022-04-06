close all;
clear;
clc;
% plot settings
FontSize = 13;
LineWidth = 1.5;
MarkerSize = 5;

surf(peaks(100))
xlabel('x(m)')
ylabel('y(m)')
zlabel('Depth(m)')
set(gca,'fontsize',FontSize)
set(gca,'fontsize',FontSize)
set(gca,'FontName','Arial','fontsize',FontSize,'Linewidth',LineWidth,'fontweight','normal')
saveas(gcf, 'Peaks', 'png')