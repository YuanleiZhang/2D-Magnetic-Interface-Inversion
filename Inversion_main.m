%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%        Inversion progress    %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%    only invert height of interface     %%%%%%%%%%%%
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
inv_M = 15;
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
% **************** 根据已知的磁化强度信息拟合出多项式函数 ***********
x = inv_x_model';
y = zeros(1, inv_N);%初始化迭代矩阵，防止有初值的影响
for i = 1 : size(coeff_M, 2)
        for j = 1 : inv_N
            y(j) = y(j)+ coeff_M(i)* x(j)^(i-1);
        end
end
% 选择用 常量 or 多项式 磁化强度
inv_M = inv_M* ones(inv_N, 1); % constant
% inv_M = y';                   % ploynamial
m_0 = inv_z_buttom/2 * ones(inv_N, 1);
% Reference model
m_ref = 0 * ones(inv_N, 1);
% Maximum number of iterations 
maxit = 30; 
% Tolerance
global TOL;
global misfit_target;
TOL = 1e-5;
misfit_target = 0.3; % \Phi_d/ Number of observation
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
N_Wm = size(m_0, 1);
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

%% Wm2 最光滑约束
W_m = W_m_2;
tic % 计时
% [recover_model, Rms] = gauss_newton_inversion(maxit, lambda, ...
%                             observation_Delta_T, x_observation, z_observation,...
%                             inv_x_left, inv_x_right, m_0, inv_z_buttom, M, Is, m_ref);

[recover_model, Rms] = gauss_newton_inversion_cool(maxit, max_lambda, num_lambda, cooling_rate,...
                                        observation_Delta_T, x_observation, ...
                                        z_observation, inv_x_left, inv_x_right,...
                                        m_0, inv_z_buttom, inv_M, Is,m_ref);                        
% for i = 1: maxit
%      [J] = compute_jacobi(x_observation, z_observation, inv_x_left, inv_x_right, m_0, inv_z_buttom, M, Is);
%      [Hax_m0, Za_m0, delta_T_m0] = magnetic_forward_2D_Guan(x_observation, z_observation, inv_x_left, inv_x_right, m_0, inv_z_buttom, M, Is);
%      % 选择合适的正则化因子
%      numsteps = 10; % 计算10个正则化因子的值
%      lambda_max = 1;
%      lambda_min = 1e0 * lambda_max;
%      lambda = logspace(log10(lambda_min),log10(lambda_max),numsteps);
%      m_try = m_0;%+
%      for j = 1:numsteps
%          g_k = - J' * (W_d' * W_d) * (observation_Delta_T - delta_T_m0) + lambda(j) * (W_m' * W_m) * (m_0 - m_ref);
%          H_k = J' * (W_d' * W_d) * J + lambda(j) * (W_m' * W_m);
%          delta_m(:,j) = pinv(H_k) * g_k;
%          m_try = m_try - delta_m(:,j);
% %          if (max(m_try) > 0 || min(m_try) < inv_z_buttom) % 上下限约束
% %              m_try = (max(m_try) + 1 - m_try)./ (max(m_try) - min(m_try)) * (inv_z_buttom);
% %          end
%          [Hax_pre(:,j), Za_pre(:,j),d_pre(:,j)] = magnetic_forward_2D_Guan(x_observation, z_observation, inv_x_left, inv_x_right, m_try, inv_z_buttom, M, Is);
%          misfit(j) = ((observation_Delta_T - d_pre(:,j))'*W_d'*W_d*(observation_Delta_T - d_pre(:,j)))/length(observation_Delta_T); %sqrt((rho_noise - d_pre(:,j))'*(rho_noise - d_pre(:,j))/(length(rho_noise)));    
%      end
%      good_misfit = min(misfit);
%      [good_misfit_row,good_misfit_column] = find(misfit == good_misfit); % 给出最小值 good_misfit 在矩阵(向量)misfit中的行号row和列号column
%      good_delta_m = delta_m(:,good_misfit_column);
%      Rms(i) = good_misfit;
%      if Rms(i) < misfit_target
%          break;
%      else
%          m_0 = m_0 - good_delta_m;
% %          if (max(m_0) > 0 || min(m_0) < inv_z_buttom) % 上下限约束
% %              m_0 = (max(m_0) + 1 - m_0)./ (max(m_0) - min(m_0)) * (inv_z_buttom);
% %          end
%      end
% end
toc
[recover_model_Hax, recover_model_Za, recover_model_delta_T]=... 
                magnetic_forward_2D_Guan(x_observation, z_observation, inv_x_left,...
                                    inv_x_right, recover_model, inv_z_buttom, inv_M, Is);

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
plot(inv_x_model, inv_M, 'r--', 'Linewidth', LineWidth)
ylabel('Magnetization(A/m)')
ax = gca;
ax.YColor = 'r';
set(gca,'fontsize',FontSize)
set(gca,'FontName','Arial','fontsize',FontSize,'Linewidth',LineWidth,'fontweight','normal')
saveas(gcf, 'Inversion result-height-cooling', 'png')

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
saveas(gcf, 'Inversion result-Fitting-cooling', 'png')