% ********************************************************* %
% ***************   Gauss_Newton Inversion  *************** %
% ***************      Cooling  Lambda      *************** %
% *******    variable :   z_up [MM*1] and [MM*1]     ****** %
% ********************************************************* %

function [recover_model, recover_M, Rms] = invert_z_up_and_M(maxit, max_lambda, num_lambda, cooling_rate,...
                                        observation_Delta_T, x_observation, ...
                                        z_observation, inv_x_left, inv_x_right,...
                                        inv_z_up_0, inv_z_buttom, M_0, Is, m_ref, M_ref, M_min, M_max)
    global W_d;
    global W_m;
    global TOL;
    global misfit_target;
    lambda = max_lambda;
    I = ones(size(inv_z_up_0, 1));  % A [MM * MM] unit matrix
    [Hax_m0, Za_m0, delta_T_m0] = magnetic_forward_2D_Guan(x_observation, z_observation,...
                                                    inv_x_left, inv_x_right, inv_z_up_0, inv_z_buttom, M_0, Is);
        Rms(1, 1) = sqrt(((observation_Delta_T - delta_T_m0)'...
                            *(W_d' * W_d)*(observation_Delta_T - delta_T_m0))/(length(observation_Delta_T)));
        disp(['Lambda = ', num2str(lambda), '       ', 'misfit = ', num2str(Rms(1, 1))]);
        misfit_last = Rms(1,1);
    for k = 1 : num_lambda
        for j = 1 : maxit
            misfit_last_last = misfit_last;
            [Hax_m0, Za_m0, delta_T_m0] = magnetic_forward_2D_Guan(x_observation, z_observation,...
                                                    inv_x_left, inv_x_right, inv_z_up_0, inv_z_buttom, M_0, Is);
            
            % Calculate Jacobian matrix  A[k] of the froward function with respect to the top verctor $model_z_up$.
            [A_k] = compute_jacobi(x_observation, z_observation,...
                                     inv_x_left, inv_x_right, inv_z_up_0, inv_z_buttom, M_0, Is);
            J_k = - 2 * A_k' * (W_d' * W_d) * (observation_Delta_T - delta_T_m0)...
                  + 2 * lambda * (W_m' * W_m) * (inv_z_up_0 - m_ref);                          % equation (8)
            H_k = 2 * A_k' * (W_d' * W_d) * A_k + 2 * lambda * (W_m' * W_m);                   % equation (11)
            delta_z_up = - pinv(H_k) * J_k;                                                        %equation (17)
            
            % Calculate Jacobian matrix  A[k] of the froward function with respect to the top verctor $ Ms $
            [B_k] = compute_jacobi_M(x_observation, z_observation,...
                                     inv_x_left, inv_x_right, inv_z_up_0, inv_z_buttom, M_0, Is);
            J_k_M = - 2 * B_k' * (W_d' * W_d) * (observation_Delta_T - delta_T_m0)...
                  + 2 * lambda * (W_m' * W_m) * (M_0 - M_ref);
            H_k_M = 2 * B_k' * (W_d' * W_d) * B_k + 2 * lambda * (W_m' * W_m);
            delta_M = - pinv(H_k_M) * J_k_M;

            inv_z_up_0 = inv_z_up_0 + delta_z_up;
            if (max(inv_z_up_0) > 0 || min(inv_z_up_0) < inv_z_buttom) % 上下限约束 -  约束基岩的上界面只能在 [0, model_z_buttom]
                    inv_z_up_0 = (max(inv_z_up_0) + 1 - inv_z_up_0)./ (max(inv_z_up_0) - min(inv_z_up_0)) * (inv_z_buttom);
            end
            M_0 = M_0 + delta_M;
            % 加入 上下限约束 - 约束岩石的物性在先验资料的磁化率范围内  -  使用归一化方法
            if (max(M_0) > M_max || min(M_0) < M_min) % 上下限约束 -  约束基岩的磁化强度只能在[M_min, M_max]
                    M_0 = M_min + ((M_0) - min(M_0))./ (max(M_0) - min(M_0)) * (M_max - M_min);
            end
            
            [Hax_m0, Za_m0, delta_T_m0] = magnetic_forward_2D_Guan(x_observation, z_observation,...
                                                    inv_x_left, inv_x_right, inv_z_up_0, inv_z_buttom, M_0, Is);
            Rms(k, j + 1) = sqrt(((observation_Delta_T - delta_T_m0)'...
                             *W_d' * W_d*(observation_Delta_T - delta_T_m0))/(length(observation_Delta_T)));
            misfit_last = Rms(k, j + 1);
            disp(['Lambda = ', num2str(lambda), '       ', 'misfit = ', num2str(Rms(k, j + 1))]);
            if misfit_last < misfit_target
                    disp('Stop, because the target misift has been achieved!');
                    break;
            elseif misfit_last > misfit_last_last
                disp('Stop, because the misift start increase!');
    %                 m_0 = m_0 + delta_m;
                lambda = cooling_rate * lambda;
                break;
            elseif abs(misfit_last - misfit_last_last) < TOL
                disp('Stop, because the misift stagnated!');
                lambda = cooling_rate * lambda;
                break;
            end
        end
        if misfit_last < misfit_target
            break;
        end
        if k == num_lambda
            disp('Stop, because the max iterative number has been achieved!');
        break;
        end

    end
    recover_M = M_0;
    recover_model = inv_z_up_0;
end
