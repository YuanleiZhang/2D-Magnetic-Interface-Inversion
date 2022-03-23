%*********************************************************************%
%******************     Calculate Jacoobi Matrix     *****************%
%*********************************************************************%
% Input parameter :
%   x - Coordinate of observation sites： ob_x  [NN*1]; Unit(m)
%   z - Coordinate of observation sites： ob_z  [NN*1]; Unit(m)
%   the left boundary of the rectangle ： x_min  [MM*1} :Unit(m)
%   the right boundary of the rectangle:  x_max  [MM*1] :Unit(m)
%   the top boundary of the rectangle  ： z_up   [MM*1] :Unit(m) 
%   the buttom boundary of the rectangl ：z_down [1*1] :Unit(m)   known number
%   magnetization ampitude ：M  [MM*1]; Unit (A/m) This is the inversion parameter
%                                        in this program to compute Jacoobi Matrix
%   magnetized direction ：I_s  [1*1]; Unit(。)

% Output : Jacoobi Matrix --->[NN*MM] 

% Author： Yuanlei Zhang
% Time :2022/03/12
function [Jacobi] = compute_jacobi_M(ob_x, ob_z, x_min, x_max, z_up, z_down, Ms, I_s)
    NN = size(ob_x, 1); MM = size(z_up, 1);
    step = 1e-5;
    Jacobi = zeros(NN, MM);
    for k = 1:MM
        m_delta = Ms;
        delta = step * abs(Ms(k));
        m_delta(k) = m_delta(k) + delta;
        Ms_delta = m_delta;
        [Hax_1, Za_1, Delta_T_1] = magnetic_forward_2D_Guan(ob_x, ob_z, x_min, x_max, z_up, z_down, Ms, I_s);
        [Hax_2, Za_2, Delta_T_2] = magnetic_forward_2D_Guan(ob_x, ob_z, x_min, x_max, z_up, z_down, Ms_delta, I_s);
        Jacobi(:,k) = (1.0/delta)*(Delta_T_2 - Delta_T_1);
    end
end