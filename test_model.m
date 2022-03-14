%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%  test_forward  %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%        test_singularity      %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author： Yuanlei Zhang
% Time :2022/03/11
close all; clear; clc

%% model parameter
x1 = -5;
x2 = 5;
z1 = 0;
z2 = -5;
M = 10;
Is = 90;
I0 = 90;
% x1 = 400;
% x2 = 600;
% z1 = -150;
% z2 = -300;
% M = 100;
% Is = 45;
% I0 = 45;

%% observation sites
site_x = -30:1:30;
% site_x = 0 : 10 : 1000;
site_z = 1;

%% calculate magnetic components ----analytical solutions
for i = 1:size(site_x, 2)
    [Hax(i,1), Za(i, 1), Delta_T(i,1)] = magnetic_forward_2D(site_x(i), site_z, x1, x2, z1, z2, M, Is);
    % Delta_T(i,1) = Hax(i,1) * cos(pi * I0 / 180) + Za(i, 1)*sin(pi * I0 / 180);
    % Delta_T(i,1) = sqrt(Hax(i,1)^2 + Za(i,1)^2);
end

for i = 1:size(site_x, 2)
    [Hax_Guan(i,1), Za_Guan(i, 1), Delta_T_Guan(i,1)] = magnetic_forward_2D_Guan(site_x(i), site_z, x1, x2, z1, z2, M, Is);
    % Delta_T(i,1) = Hax(i,1) * cos(pi * I0 / 180) + Za(i, 1)*sin(pi * I0 / 180);
    % Delta_T(i,1) = sqrt(Hax(i,1)^2 + Za(i,1)^2);
end
    [Hax_vector, Za_vector, Delta_T_vector] = magnetic_forward_2D(site_x, site_z, x1, x2, z1, z2, M, Is);
% for i = 1:size(site_x, 2)
%     f_Hax = @(x,z) (site_x(i) - x) ./ (((site_x(i) - x).*(site_x(i) - x) + (site_z - z).*(site_z - z)).^(3/2));
%     I_Hax(i, 1) = -1e2 * M * sin(pi * Is / 180) * integral2(f_Hax, x1, x2, z1, z2);
%     f_Za = @(x,z) (site_z - z) ./ (((site_x(i) - x).*(site_x(i) - x) + (site_z - z).*(site_z - z)).^(3/2));
%     I_Za(i, 1) = -1e2 * M * sin(pi * Is / 180) * integral2(f_Za, x1, x2, z1, z2);
% end

%% plot 
% plot settings
FontSize = 15;
LineWidth = 1.5;
MarkerSize = 5;

figure('Position',[200,100,1250,600])
subplot(2,2,1)
plot(site_x, zeros(size(site_x,2)), 'k-','LineWidth', LineWidth)
hold on
plot(site_x, ones(size(site_x,2)), 'r-.','LineWidth', LineWidth)
xlim([site_x(1), site_x(end)])
ylim([z2 + 0.1*z2, site_z + abs(0.1*z2)])
rectangle('position',[x1, z2, x2 - x1, abs(z2 - z1)], ...
          'EdgeColor', 'k', ...
          'LineWidth', LineWidth);
quiver(x1 + 0.1*(x2 - x1), z1 + 0.1*(z2 - z1),...
        0.7*sqrt((z2 - z1)*(z2 - z1) + (x2 - x1)*(x2 - x1))* cos(pi * Is / 180),...
        -0.7*sqrt((z2-z1)*(z2-z1) + (x2-x1)*(x2-x1))* sin(pi * Is / 180),...
        'MaxHeadSize', 23, 'color', 'b', 'LineWidth',LineWidth)
set(gca,'xaxislocation','top','yaxislocation','left');
xlabel('x(m)')
ylabel('Depth(m)')
set(gca,'fontsize',FontSize)
set(gca,'FontName','Arial','fontsize',FontSize,'Linewidth',LineWidth,'fontweight','normal')      
subplot(2,2,2)
plot(site_x, Hax, 'k-','LineWidth', LineWidth)
hold on
plot(site_x, Hax_Guan, 'ro','LineWidth', LineWidth)
plot(site_x, Hax_vector, 'b^','LineWidth', LineWidth)
xlabel('x(m)')
ylabel('H_{ax}(nT)')
legendon = legend('Analytic solution of Liu(2020)','Analytic solution of Guan(2005)','Location','NorthWest');
set(legendon,'box','off')
set(gca,'fontsize',FontSize)
set(gca,'FontName','Arial','fontsize',FontSize,'Linewidth',LineWidth,'fontweight','normal')
subplot(2,2,3)
plot(site_x, Za, 'k-','LineWidth', LineWidth)
hold on
plot(site_x, Za_Guan, 'or','LineWidth', LineWidth)
plot(site_x, Za_vector, 'b^','LineWidth', LineWidth)
xlabel('x(m)')
ylabel('Z_{a}(nT)')
set(gca,'fontsize',FontSize)
set(gca,'FontName','Arial','fontsize',FontSize,'Linewidth',LineWidth,'fontweight','normal')
subplot(2,2,4)
plot(site_x, Delta_T, 'k-','LineWidth', LineWidth)
hold on
plot(site_x, Delta_T_Guan, 'or','LineWidth', LineWidth)
xlabel('x(m)')
ylabel('\Delta T(nT)')
set(gca,'fontsize',FontSize)
set(gca,'FontName','Arial','fontsize',FontSize,'Linewidth',LineWidth,'fontweight','normal')

saveas(gcf, 'test', 'png')