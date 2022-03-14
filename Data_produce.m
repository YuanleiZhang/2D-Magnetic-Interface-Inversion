%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%  Data produce  %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author： Yuanlei Zhang
% Time :2022/03/12
close all;
clear;
clc;

%% Data produce  
% magnetization amplitude; Unit : (A/m)
M = 10;
% Magnetized direction; Unit : degree(。)
Is = 90;
I0 = 90;
% range of horizontal direction;  Unit : (m)
x_start = 0;
x_end = 100;
% average depth of interface 
average_depth = -20;
% Divided numbers : N
N = 100;
diff = (x_end - x_start) / N;
x_model = (x_start + diff/2 : diff : x_end - diff/2)';
x_model_left = (x_start : diff : x_end - diff)';
x_model_right = (x_start + diff : diff : x_end)';
% Random data in peaks(N) : n
n = 36; 
grid = peaks(N);
% Original model 
model_z_up = average_depth - 3 * grid(:, n) + 5;
% Buttom depth of interface
model_z_buttom = floor(min(model_z_up*0.1))*10;

% Observation sites : x 
x_observation = (x_start : 1 : x_end)';

% Observation height
z_observation = 1;

% Caculate magntic field Hax, Za and Delta_T
for i = 1: size(x_observation, 1)
    Hax(i, 1) = 0;
    Za(i, 1) = 0;
    Delta_T(i, 1) = 0;
    for j = 1 : size(model_z_up, 1)
        [Hax_temp, Za_temp, Delta_T_temp] = magnetic_forward_2D(x_observation(i), z_observation, ...
                              x_start + diff*(j - 1), x_start + diff*(j), model_z_up(j),...
                              model_z_buttom, M, Is);
        Hax(i,1) = Hax(i,1) + Hax_temp;
        Za(i,1) = Za(i,1) + Za_temp;
%         Delta_T_temp = Hax_temp * cos(pi * I0 / 180) + Za_temp*sin(pi * I0 / 180);
        Delta_T(i, 1) = Delta_T(i, 1) + Delta_T_temp;
    end
end
% Calculate magnetic field use vector directly
[Hax_vector, Za_vector,Delta_vector]  = magnetic_forward_2D(x_observation, z_observation, ...
                               x_model_left, x_model_right, model_z_up,...
                               model_z_buttom, M, Is);
%*********************************************************************************************************************
% Hax_vector = zeros(size(x_observation, 1), 1);
% Za_vector = zeros(size(x_observation, 2), 1);
% Delta_T_vector = zeros(size(x_observation, 1), 1);
% for i = 1:size(model_z_up, 1)
%     [Hax_vector_temp, Za_vector_temp, Delta_T_vector_temp] = magnetic_froward_2D(x_observation, z_observation, ...
%                               x_start + diff*(i - 1), x_start + diff*(i), model_z_up(i),...
%                               model_z_buttom, M, Is);
%     Hax_vector = Hax_vector + Hax_vector_temp;
%     Za_vector = Za_vector + Za_vector_temp;
%     Delta_T_vector = Delta_T_vector + Delta_T_vector_temp; 
% end
%*********************************************************************************************************************

% Add Noise into Delta_T
epsilon = 0.05;
noise = 2*(0.5-rand(size(Delta_T)))* epsilon .* Delta_T;
observation_Delta_T = Delta_T + noise;
save('magnetic_responce','x_model', 'x_model_left', 'x_model_right', 'model_z_up','model_z_buttom','x_observation', 'z_observation', 'Hax', 'Za', 'Delta_T', 'observation_Delta_T');



%% plot section
% plot settings
FontSize = 13;
LineWidth = 1.5;
MarkerSize = 5;

% plot_model
figure('Position',[100,200,1200,300])
subplot(1,2,1)
plot(x_model, model_z_up, 'k-', 'LineWidth', LineWidth)
xlabel('x(m)')
ylabel('Depth(m)')
title('Orignal model')
set(gca,'fontsize',FontSize)
set(gca,'FontName','Arial','fontsize',FontSize,'Linewidth',LineWidth,'fontweight','normal')
subplot(1,2,2)
bar(x_model, model_z_up - model_z_buttom)
xlabel('x(m)')
ylabel('Depth(m)')
title('Divided the Orignal Model into Quadrilaterals')
set(gca,'fontsize',FontSize)
set(gca,'FontName','Arial','fontsize',FontSize,'Linewidth',LineWidth,'fontweight','normal')
ylim([0, -model_z_buttom])
for i = 1 : abs(model_z_buttom / 10) + 1
    Ylabel(i, 1) = (model_z_buttom + (i - 1)*10);
end
set(gca,'YTickLabel', num2str(Ylabel))
saveas(gcf, 'Orignal model', 'png')
%%%%%%%%%%%%%%%%%%%%%%%% plot forward responce   %%%%%%%%%%%%%%%%%%%%
figure('Position',[200,100,750,500]) 
subplot(2,2,1)
plot(x_observation, Hax, 'k-','LineWidth', LineWidth)
hold on
plot(x_observation, Hax_vector, 'ro','LineWidth', LineWidth)
xlabel('x(m)')
ylabel('H_{ax}(nT)')
set(gca,'fontsize',FontSize)
set(gca,'FontName','Arial','fontsize',FontSize,'Linewidth',LineWidth,'fontweight','normal')
subplot(2,2,2)
plot(x_observation, Za, 'k-','LineWidth', LineWidth)
hold on
plot(x_observation, Za_vector, 'ro','LineWidth', LineWidth)
xlabel('x(m)')
ylabel('Z_{a}(nT)')
set(gca,'fontsize',FontSize)
set(gca,'FontName','Arial','fontsize',FontSize,'Linewidth',LineWidth,'fontweight','normal')
subplot(2,2,3)
plot(x_observation, Delta_T, 'k-','LineWidth', LineWidth)
xlabel('x(m)')
ylabel('\Delta T(nT)')
set(gca,'fontsize',FontSize)
set(gca,'FontName','Arial','fontsize',FontSize,'Linewidth',LineWidth,'fontweight','normal')
subplot(2,2,4)
plot(x_observation, observation_Delta_T, 'k-','LineWidth', LineWidth)
xlabel('x(m)')
ylabel(['\Delta T within ', num2str(epsilon*100), '% noise (nT)'])
set(gca,'fontsize',FontSize)
set(gca,'FontName','Arial','fontsize',FontSize,'Linewidth',LineWidth,'fontweight','normal')
saveas(gcf, 'Observation data', 'png')