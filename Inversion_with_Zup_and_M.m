%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%        Inversion progress    %%%%%%%%%%%%%%%%%%%%%%%%
%%%%   invert the height of interface and magnetization constrast    %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Reference : Marlon C. Hidalgo-Gato,et al. Magnetic amplitude inversion
%%%%             for depth-to-basement and apparent magnetization-intensity
%%%%             estimates, Geophysics,(2021),86(1):J1-J11,10.1190/GEO2019-0726.1

% Author： Yuanlei Zhang
% Time :2022/03/11

close all;
clear;
clc;
load magnetic_responce.mat

%% %%%%%%%%%%%%%%%%%   Inversion progress    %%%%%%%%%%%%%%%%% 
% Inversion parameter
% Magnetized direction; Unit : degree(。)
Is = 90;
I0 = 90;
% range of inversion horizontal direction;  Unit : (m)
x_start = 0;
x_end = 100;
% Number of inversion sections
inv_N = 150;
inv_diff = (x_end - x_start)/inv_N;
inv_x_model = (x_start + inv_diff/2 : inv_diff :  x_end - inv_diff/2)';
inv_x_left = (x_start : inv_diff : x_end - inv_diff)';
inv_x_right = (x_start + inv_diff : inv_diff : x_end)';
% Intital model (the buttom of Rock is known) equal to z_buttom
inv_z_buttom = model_z_buttom;
% Reference model
m_ref = 0 * ones(inv_N, 1);
M_ref = 0 * ones(inv_N, 1);
% Maximum number of iterations 
maxit = 20; 
% Tolerance
global TOL;
global misfit_target;
TOL = 1e-3;
misfit_target = 1; % \Phi_d/ Number of observation
% Constant lambda for 'gauss_newton_inversion.m'
lambda = 1;
% Intital number of regularization factor 
max_lambda = 100;
num_lambda = 15;
% Cooling coefficient 
cooling_rate = 0.6;
% Calculate data covariance matrix -Wd- and model covariance matrix -Wm- 
% Wd
N_Wd = size(observation_Delta_T, 1);
N_Wm = inv_N;
global W_d;
global W_m;
W_d = zeros(N_Wd, N_Wd);
% relative_error
relative_error = 0.05;
for i = 1:N_Wd
   W_d(i,i) = 1/(relative_error * abs(observation_Delta_T(i)));
end
% W_m_1
W_m_1 = zeros(N_Wm, N_Wm);
W_m_1(1,1) = 0;
W_m_1(1,2) = 0;
for i = 2:N_Wm
    W_m_1(i,i) = 1;
    W_m_1(i,i-1) = -1;
end
% W_m_2
W_m_2 = zeros(N_Wm, N_Wm);
for i = 3:N_Wm
    W_m_2(i,i) = 1;
    W_m_2(i,i-1) = -2;
    W_m_2(i,i-2) = 1;
end

%% Start inversion
% inversion model magnetization M --> M and model_z_up --> m_0
W_m = W_m_1;
M_0 = 1 * ones(inv_N, 1);
m_0 = inv_z_buttom/2 * ones(inv_N, 1);
M_ref_min = 5;
M_ref_max = 20;
tic % 计时
% [recover_model, recover_M, Rms] = invert_z_up_and_constant_M(maxit, max_lambda, num_lambda, cooling_rate,...
%                                         observation_Delta_T, x_observation, ...
%                                         z_observation, inv_x_left, inv_x_right,...
%                                         m_0, inv_z_buttom, M_0, Is, m_ref);

[recover_model, recover_M, Rms] = invert_z_up_and_M(maxit, max_lambda, num_lambda, cooling_rate,...
                                        observation_Delta_T, x_observation, ...
                                        z_observation, inv_x_left, inv_x_right,...
                                        m_0, inv_z_buttom, M_0, Is, m_ref, M_ref, M_ref_min, M_ref_max);
                                    
toc
[recover_model_Hax, recover_model_Za, recover_model_delta_T]=... 
                magnetic_forward_2D_Guan(x_observation, z_observation, inv_x_left,...
                            inv_x_right, recover_model, inv_z_buttom, recover_M, Is);
                                
%% plot section
% plot settings
FontSize = 13;
LineWidth = 1.5;
MarkerSize = 5;

%%%%%%%%%%%%%%%%%%%%%%%    Inversion model    %%%%%%%%%%%%%%%%%%%%%%%%
figure('Position',[200,400,1200,300])
subplot(1,2,1)
yyaxis left;
plot(x_model, model_z_up, 'k-', 'Linewidth', LineWidth)
xlabel('x(m)')
ylabel('Depth(m)')
title('Orignal model')
ax = gca;
ax.YColor = 'k';
yyaxis right;
plot(x_model, M, 'r--', 'Linewidth', LineWidth)
hold on
plot(x0, y0, 'bd', 'Linewidth', LineWidth)
ylabel('Magnetization(A/m)')
ax = gca;
ax.YColor = 'r';
set(gca,'fontsize',FontSize)
set(gca,'FontName','Arial','fontsize',FontSize,'Linewidth',LineWidth,'fontweight','normal')
subplot(1,2,2)
yyaxis left;
plot(inv_x_model, recover_model, 'k-', 'Linewidth', LineWidth)
xlabel('x(m)')
ylabel('Depth(m)')
title('Orignal model')
ax = gca;
ax.YColor = 'k';
yyaxis right;
plot(inv_x_model, recover_M, 'r--', 'Linewidth', LineWidth)
ylabel('Magnetization(A/m)')
ax = gca;
ax.YColor = 'r';
set(gca,'fontsize',FontSize)
set(gca,'FontName','Arial','fontsize',FontSize,'Linewidth',LineWidth,'fontweight','normal')
saveas(gcf, 'Inversion result-delta(M and z)-height-cooling', 'png')

 %****************  Plot data fitting and misfit   **********************
figure('Position',[300,300,950,300]) 
subplot(1,2,1)
plot(x_observation, observation_Delta_T, 'k^-','LineWidth', LineWidth)
hold on
plot(x_observation, recover_model_delta_T, 'r>-','LineWidth', LineWidth)
xlabel('x(m)')
ylabel('\Delta T (nT)')
legendon = legend('Observation Data', '\Delta_T of Recover model', 'Location', 'best');
set(legendon,'box','off')
set(gca,'fontsize',FontSize)
set(gca,'FontName','Arial','fontsize',FontSize,'Linewidth',LineWidth,'fontweight','normal')
subplot(1,2,2)
Rms_temp = Rms';
Rms_temp(find(Rms_temp==0))=[];
plot(Rms_temp,'k+-','LineWidth', LineWidth)
% RMS2_end = Rms2(end);
xlabel('Iteration (Times)')
ylabel('RMS(%)')
set(gca,'fontsize',FontSize)
set(gca,'fontsize',FontSize)
set(gca,'FontName','Arial','fontsize',FontSize,'Linewidth',LineWidth,'fontweight','normal')
saveas(gcf, 'Inversion result-delta(M and z)-Fitting-cooling', 'png')