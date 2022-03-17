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
% Maximum number of iterations 
maxit = 20; 
% Tolerance
global TOL;
global misfit_target;
TOL = 1e-5;
misfit_target = 0.3; % \Phi_d/ Number of observation
% Constant lambda for 'gauss_newton_inversion.m'
lambda = 1;
mu = 0.2;
% Intital number of regularization factor 
max_lambda = 100;
num_lambda = 20;
% Cooling coefficient 
cooling_rate = 0.6;
% Calculate data covariance matrix -Wd- and model covariance matrix -Wm- 
% Wd
N_Wd = size(observation_Delta_T, 1);
N_Wm = size(inv_N, 1);
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
W_m = W_m_2;
M_0 = 1;
m_0 = inv_z_buttom/2 * ones(inv_N, 1);

tic % 计时
lambda = max_lambda;
I = ones(size(m_0, 1));  % A [MM * MM] unit matrix
[Hax_m0, Za_m0, delta_T_m0] = magnetic_forward_2D_Guan(x_observation, z_observation,...
                                                inv_x_left, inv_x_right, m_0, inv_z_buttom, M_0, Is);
    Rms(1, 1) = ((observation_Delta_T - delta_T_m0)'...
                        *(W_d' * W_d)*(observation_Delta_T - delta_T_m0))/(length(observation_Delta_T));
    disp(['Lambda = ', num2str(lambda), '       ', 'misfit = ', num2str(Rms(1, 1))]);
for k = 1 : num_lambda
    for j = 1 : maxit
        [Hax_m0_k, Za_m0_k, delta_T_m0_k] = magnetic_forward_2D_Guan(x_observation, z_observation,...
                                                inv_x_left, inv_x_right, m_0, inv_z_buttom, 1, Is);
        %  update of m_{k+1}   equation(16)
        M_next = (observation_Delta_T' * delta_T_m0_k) / (delta_T_m0_k' * delta_T_m0_k);
        % Calculate Jacobian matrix  A[k] of the froward function with respect to the top verctor model_z_up.
        [A_k] = compute_jacobi(x_observation, z_observation,...
                                 inv_x_left, inv_x_right, m_0, inv_z_buttom, M_0, Is);
        J_k = -2 * A_k' * (W_d' * W_d) * (observation_Delta_T - M_0 * delta_T_m0) + 2 * lambda * (W_m' * W_m) * (m_0 - m_ref);  % equation (8)
        H_k = 2 * A_k' * (W_d' * W_d) * A_k + 2 * lambda * (W_m' * W_m);   % equation (11)
        delta_m = - pinv(H_k) * J_k;   %equation (17)
        if (max(delta_m) > 0 || min(delta_m) < inv_z_buttom) % 上下限约束
                delta_m = (max(delta_m) + 1 - delta_m)./ (max(delta_m) - min(delta_m)) * (inv_z_buttom);
        end
        M_0 = M_next;
        m_0 = m_0 - delta_m;
        if (max(m_0) > 0 || min(m_0) < inv_z_buttom) % 上下限约束
                m_0 = (max(m_0) + 1 - m_0)./ (max(m_0) - min(m_0)) * (inv_z_buttom);
        end
        [Hax_m0, Za_m0, delta_T_m0] = magnetic_forward_2D_Guan(x_observation, z_observation,...
                                                inv_x_left, inv_x_right, m_0, inv_z_buttom, M_0, Is);
        Rms(k, j + 1) = ((observation_Delta_T - delta_T_m0)'...
                         *(W_d' * W_d)*(observation_Delta_T - delta_T_m0))/(length(observation_Delta_T));
        misfit_last = Rms(k, j + 1); 
        disp(['Lambda = ', num2str(lambda), '       ', 'misfit = ', num2str(Rms(k, j + 1))]);
        if Rms(k, j + 1) < misfit_target
                disp('Stop, because the target misift has been achieved!');
                break;
        elseif Rms(k, j + 1) > Rms(k, j)
            disp('Stop, because the misift start increase!');
%                 m_0 = m_0 + delta_m;
            lambda = cooling_rate * lambda;
            break;
        elseif abs(Rms(k, j + 1) - Rms(k, j)) < TOL
            disp('Stop, because the misift stagnated!');
            lambda = cooling_rate * lambda;
            break;
        end
        j = j+1;
    end
    if misfit_last < misfit_target
            break;
    end
    if k == num_lambda
        disp('Stop, because the max iterative number has been achieved!');
        break;
    end
end
recover_M = M_0*ones(inv_N,1);
recover_model = m_0;
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