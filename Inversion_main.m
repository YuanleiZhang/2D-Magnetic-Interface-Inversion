%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%        Inversion progress    %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%  only inversion height of interface     %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author： Yuanlei Zhang
% Time :2022/03/11
close all;
clear;
clc;
load magnetic_responce.mat
%% %%%%%%%%%%%%%%%%%   Inversion progress    %%%%%%%%%%%%%%%%% 
% Inversion parameter
% magnetization amplitude; Unit : (A/m)
M = 10;
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

m_0 = model_z_buttom/3 * ones(inv_N, 1);
% Reference model
m_ref = model_z_buttom/2 * ones(inv_N, 1);
% Maximum number of iterations 
maxit = 20; 
% Tolerance
TOL = 1e-3; 
% Intital number of regularization factor 
lambda = 1e-1;
% Cooling coefficient 
c = 0.6;

% Calculate data covariance matrix -Wd- and model covariance matrix -Wm- 
% Wd
N_Wd = size(observation_Delta_T, 1);
N_Wm = size(m_0, 1);
W_d = zeros(N_Wd, N_Wd);
% relative_error
relative_error = 0.02;
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

%% Wm2 最光滑约束
W_m = W_m_2;
tic % 计时
for i = 1:maxit 
    [J] = compute_jacobi(x_observation, z_observation, inv_x_left, inv_x_right, m_0, inv_z_buttom, M, Is);
    g_k = - J' * (W_d' * W_d) * (observation_Delta_T - magnetic_froward_2D(x_observation, z_observation,...
                                         inv_x_left, inv_x_right, m_0, inv_z_buttom, M, Is)) + lambda * (W_m' * W_m) * (m_0 - m_ref);
    H_k = J' * (W_d' * W_d) * J + lambda * (W_m' * W_m);
    delta_m = pinv(H_k) * g_k;
    [Hax_pre, Za_pre,d_pre] = magnetic_froward_2D(x_observation, z_observation, inv_x_left, inv_x_right, m_0, inv_z_buttom, M, Is);
    Rms(i) = norm(W_d*(observation_Delta_T - d_pre), 2)/(length(observation_Delta_T));
    if Rms < 1  
        break;
    else
        m_0 = m_0 - delta_m;
    end
end
% for i = 1: maxit 
%      [J] = compute_jacobi(x_observation, z_observation, inv_x_left, inv_x_right, m_0, inv_z_buttom, M, Is);
%      % 选择合适的正则化因子
%      numsteps = 10; % 计算10个正则化因子的值
%      lambda_max = 100;
%      lambda_min = 1e-4 * lambda_max;
%      lambda = logspace(log10(lambda_min),log10(lambda_max),numsteps);
%      m_try = m_0;
%      for j = 1:numsteps
%           g_k = - J' * (W_d' * W_d) * (observation_Delta_T - magnetic_froward_2D(x_observation, z_observation,...
%                                         inv_x_left, inv_x_right, m_0, inv_z_buttom, M, Is))...
%                 + lambda(j) * (W_m' * W_m) * (m_0 - m_ref);
%           H_k = J' * (W_d' * W_d) * J + lambda(j) * (W_m' * W_m);
%          delta_m(:,j) = pinv(H_k) * g_k;
%          m_try = m_try - delta_m(:,j);
%          [Hax_pre(:,j), Za_pre(:,j),d_pre(:,j)] = magnetic_froward_2D(x_observation, z_observation, inv_x_left, inv_x_right, m_try, inv_z_buttom, M, Is);
%          misfit(j) = norm(observation_Delta_T - d_pre(:,j))/norm(observation_Delta_T); %sqrt((rho_noise - d_pre(:,j))'*(rho_noise - d_pre(:,j))/(length(rho_noise)));    
%      end
%      good_misfit = min(misfit);
%      [good_misfit_row,good_misfit_column] = find(misfit == good_misfit); % 给出最小值 good_misfit 在矩阵(向量)misfit中的行号row和列号column
%      good_delta_m = delta_m(:,good_misfit_column);
%      Rms2(i) = good_misfit; 
%      if Rms2 < TOL
%          break;
%      else
%          m_0 = m_0 - good_delta_m;
%      end
% end
toc
recover_model = m_0;
recover_model_responce = magnetic_froward_2D(x_observation, z_observation, inv_x_left, inv_x_right, recover_model, inv_z_buttom, M, Is);

%% plot section
% plot settings
FontSize = 13;
LineWidth = 1.5;
MarkerSize = 5;

%%%%%%%%%%%%%%%%%%%%%%%    Inversion model    %%%%%%%%%%%%%%%%%%%%%%%%
figure('Position',[200,400,1200,300])
subplot(1,2,1)
bar(x_model, model_z_up - model_z_buttom)
xlabel('x(m)')
ylabel('Depth(m)')
title('Orignal Model')
set(gca,'fontsize',FontSize)
set(gca,'FontName','Arial','fontsize',FontSize,'Linewidth',LineWidth,'fontweight','normal')
ylim([0, -model_z_buttom])
for i = 1 : abs(model_z_buttom / 10) + 1
    Ylabel(i, 1) = (model_z_buttom + (i - 1)*10);
end
set(gca,'YTickLabel', num2str(Ylabel))
subplot(1,2,2)
plot(recover_model, 'k-','LineWidth', LineWidth)
xlabel('x(m)')
ylabel('Depth(m)')
title('Recover Model')
set(gca,'fontsize',FontSize)
set(gca,'FontName','Arial','fontsize',FontSize,'Linewidth',LineWidth,'fontweight','normal')
% ylim([0, -model_z_buttom])
% for i = 1 : abs(model_z_buttom / 10) + 1
%     Ylabel(i, 1) = (model_z_buttom + (i - 1)*10);
% end
% set(gca,'YTickLabel', num2str(Ylabel))

 %****************  Plot data fitting and misfit   **********************
figure('Position',[300,300,850,300]) 
subplot(1,2,1)
plot(x_observation, observation_Delta_T, 'k^-','LineWidth', LineWidth)
hold on
plot(x_observation, recover_model_responce, 'r>-','LineWidth', LineWidth)
xlabel('x(m)')
ylabel('\Delta T (nT)')
legendon = legend('Observation Data', 'Recover model magnetic responce');
set(legendon,'box','off')
set(gca,'fontsize',FontSize)
subplot(1,2,2)
plot(Rms,'k+-','LineWidth', LineWidth)
% RMS2_end = Rms2(end);
xlabel('Iteration (Times)')
ylabel('RMS(%)')
set(gca,'fontsize',FontSize)