%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%   2D magnetic forward -- Guan    %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input parameter :
%   x - Coordinate of observation sites： x  [NN*1]; Unit(m)
%   z - Coordinate of observation sites： z  [NN*1] or [1*1]; Unit(m)
%   the left boundary of the rectangle ： x_min  [MM*1} : Unit(m)
%   the right boundary of the rectangle:  x_max  [MM*1] : Unit(m)
%   the top boundary of the rectangle  ： z_up   [MM*1] : Unit(m)
%   the buttom boundary of the rectangl ：z_down [MM*1] or [1*1]: Unit(m)   known number
%   magnetization ampitude ：M  [1*1] or [MM*1]; Unit (A/m)
%   magnetized direction ：I_s  [1*1]; Unit(。)
% Output :
% Horizontal component of magntic field： Hax [NN*1]; Unit : nT
% Vertical component of magntic field  ： Za  [NN*1]; Unit : nT
% Magnetic anomaly  ： delta_T  [NN*1]; Unit : nT

% Author： Yuanlei Zhang
% Time :2022/03/11

function [Hax, Za, delta_T] = magnetic_forward_2D_Guan(x, z, x_min, x_max, z_up, z_down, M, I_s)
    C = 1e-7 * 1e9;
    N_x = size(x, 1);
    N_z = size(z, 1);
    if (N_x ~= N_z) 
        if(N_z == 1 && N_x ~= 1)
            z = z * ones(N_x, 1);
        elseif (N_x == 1 && N_z ~= 1)
            x = x * ones(N_z, 1);
        else
            disp('The dimensions of the input observations are inconsistent !')
            quit;   % stop program 
        end
    end
    N_x_min = size(x_min, 1);
    N_x_max = size(x_max, 1);
    N_z_up = size(z_up, 1);
    N_z_down = size(z_down, 1);
    
    assert(N_x_min == N_x_max);
    if (N_z_up ~= N_z_down) 
        if(N_z_down == 1 && N_z_up  ~= 1)
            z_down = z_down * ones(N_z_up, 1);
            N_z_down  = size(z_down, 1);
        elseif (N_z_up == 1 && N_z_down ~= 1)
            z_up = z_up * ones(N_z_down, 1);
        else
            disp('The dimensions of the input prism are inconsistent !')
            quit;   % stop program 
        end
    end
    assert(N_x_min == N_z_down)
    N_M = size(M, 1);
    if (N_M == 1)
        M = M * ones(N_z_up);
    elseif (N_M ~= N_z_up)
        disp('The dimensions of the input prism are inconsistent !')
    end
    Hax = zeros(N_x, 1);
    Za = zeros(N_x, 1);
    delta_T = zeros(N_x, 1);
    for  i = 1 : N_x_min
        b = 0.5 * abs(x_max(i) - x_min(i));
        l = 0.5 * abs(z_down(i) - z_up(i));
        x0 = 0.5 * (x_min(i) + x_max(i));
        z0 = 0.5 * (z_up(i) + z_down(i));
        x1 = x - x0 + b;
        x2 = x - x0 - b;
        z1 = z0 - z - l;   % h 
        z2 = z0 - z + l;   % h + 2l
        %%% analytical solution 1 --- Equation (3-1-89) refered to Zhining Guan(2005) 
        A = 0.5 * log(((z1.*z1 + x2.*x2).*(z2.*z2 + x1.*x1))./((z1.*z1 + x1.*x1).*(z2.*z2 + x2.*x2)));
        F = (atan((x1)./(z1)) - atan((x2)./(z1))) - (atan((x1)./(z2)) - atan((x2)./(z2)));
        Hax_temp = C * 2 * M(i) * (sin(pi* I_s / 180) .* A - cos(pi* I_s / 180) .* F);
        Za_temp = C * 2 * M(i) * (cos(pi * I_s / 180) .* A + sin(pi * I_s / 180) .* F);
        delta_T_temp = C * 2 * M(i) * (cos(pi * (90 - 2 * I_s) / 180) .* A - sin(pi * (90 - 2 * I_s) / 180) .* F);
        %%% analytical solution 2 --- Equation (2.1 - 2.6) refered to Suang liu(2020) 
%         E = log(((z2.*z2 + x1.*x1).* (z1.*z1 + x2.*x2))./(((z1.*z1 + x1.*x1)).*(z2.*z2 + x2.*x2)));
%         F1 = atan((2.*b.*z1)./(z1.*z1 + (x1 - b).*(x1 - b) - b.*b)) - atan((2.*b.*z2)./(z2.*z2 +(x1 - b).*(x1 - b) - b.*b));
%         Hax_temp = C .* 2 .* M .* (0.5 .* sin(pi .* I_s ./ 180) .* E - cos (pi .* I_s ./ 180).* F1);
%         Za_temp = C .* 2 .* M .* (0.5 .* cos(pi .* I_s ./ 180) .* E + sin (pi .* I_s ./ 180).* F1);
%         delta_T_temp = Hax_temp .* cos(pi .* I_s ./ 180) + Za_temp .* sin(pi .* I_s ./ 180);
        Hax = Hax + Hax_temp;
        Za = Za + Za_temp;
        delta_T = delta_T + delta_T_temp;
    end
end
