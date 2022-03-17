% ********************************************************* %
% ***************   Gauss_Newton Inversion  *************** %
% ***************      Cooling  Lambda      *************** %
% ********************************************************* %

function [recover_model, Rms] = gauss_newton_inversion_cool(maxit, max_lambda, n_lambda, cooling_rate,...
                                        observation_Delta_T, x_observation, ...
                                        z_observation, inv_x_left, inv_x_right,...
                                        m_0, inv_z_buttom, M, Is,m_ref)
    global W_d;
    global W_m;
    global TOL;
    global misfit_target;
    lambda = max_lambda;
    [Hax_m0, Za_m0, delta_T_m0] = magnetic_forward_2D_Guan(x_observation, z_observation,...
                                                inv_x_left, inv_x_right, m_0, inv_z_buttom, M, Is);
    Rms(1, 1) = ((observation_Delta_T - delta_T_m0)'* W_d'...
                        *W_d*(observation_Delta_T - delta_T_m0))/(length(observation_Delta_T));
    disp(['Lambda = ', num2str(lambda), '       ', 'misfit = ', num2str(Rms(1, 1))]);
    misfit_last = Rms(1,1);
    for k = 1 : n_lambda
        for i = 1:maxit
            misfit_last_last = misfit_last;
            [J] = compute_jacobi(x_observation, z_observation,...
                                 inv_x_left, inv_x_right, m_0, inv_z_buttom, M, Is);
            g_k = - J' * (W_d' * W_d) * (observation_Delta_T - delta_T_m0) ...
                  + lambda * (W_m' * W_m) * (m_0 - m_ref);
            H_k = J' * (W_d' * W_d) * J + lambda * (W_m' * W_m);
            delta_m = pinv(H_k) * g_k;
            m_0 = m_0 - delta_m;
            if (max(m_0) > 0 || min(m_0) < inv_z_buttom) % 上下限约束
                m_0 = (max(m_0) + 1 - m_0)./ (max(m_0) - min(m_0)) * (inv_z_buttom);
            end
            [Hax_m0, Za_m0, delta_T_m0] = magnetic_forward_2D_Guan(x_observation, z_observation,...
                                            inv_x_left, inv_x_right, m_0, inv_z_buttom, M, Is);
            Rms(k, i + 1) = ((observation_Delta_T - delta_T_m0)'*W_d'...
                            *W_d*(observation_Delta_T - delta_T_m0))/(length(observation_Delta_T));
            misfit_last = Rms(k, i + 1);
            disp(['Lambda = ' num2str(lambda) '       ' 'misfit = ' num2str(Rms(k, i + 1))]);
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
        if k == n_lambda
            disp('Stop, because the max iterative number has been achieved!');
            break;
        end
    end   
    recover_model = m_0;
end
