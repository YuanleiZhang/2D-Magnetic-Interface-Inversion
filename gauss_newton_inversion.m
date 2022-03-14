% ********************************************************* %
% ***************   Gauss_Newton Inversion  *************** %
% ***************      Constant Lambda      *************** %
% ********************************************************* %

function [recover_model, Rms] = gauss_newton_inversion(maxit, lambda, observation_Delta_T, x_observation, ...
                                                z_observation, inv_x_left, inv_x_right, m_0, inv_z_buttom, M, Is,m_ref)
    global W_d;
    global W_m;
    global TOL;
    global misfit_target;
    for i = 1:maxit
        [J] = compute_jacobi(x_observation, z_observation, inv_x_left, inv_x_right, m_0, inv_z_buttom, M, Is);
        [Hax_m0, Za_m0, delta_T_m0] = magnetic_forward_2D_Guan(x_observation, z_observation, inv_x_left, inv_x_right, m_0, inv_z_buttom, M, Is);
        g_k = - J' * (W_d' * W_d) * (observation_Delta_T - delta_T_m0) + lambda * (W_m' * W_m) * (m_0 - m_ref);
        H_k = J' * (W_d' * W_d) * J + lambda * (W_m' * W_m);
        delta_m = pinv(H_k) * g_k;
        [Hax_pre, Za_pre, d_pre] = magnetic_forward_2D_Guan(x_observation, z_observation, inv_x_left, inv_x_right, m_0 - delta_m, inv_z_buttom, M, Is);
        Rms(i) = ((observation_Delta_T - d_pre)'*W_d'*W_d*(observation_Delta_T - d_pre))/(length(observation_Delta_T));
        if Rms(i) < misfit_target
            break;
        else
            m_0 = m_0 - delta_m;
        end
    end
    recover_model = m_0;
end
