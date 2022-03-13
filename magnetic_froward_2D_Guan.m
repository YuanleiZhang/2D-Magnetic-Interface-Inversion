%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%   2D magnetic forward -- Guan    %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Input parameter: 
% 观测点： x, z
% 矩形水平边界： x_min, x_max
% 矩形垂直边界： z_up, z_down
% 总磁化强度大小：M
% 磁化倾角：I_s
% Output :
% 水平方向分量： Hax / Unit : nT
% 垂直方向分量： Za / Unit : nT
function [Hax, Za, delta_T] = magnetic_froward_2D_Guan(x, z, x_min, x_max, z_up, z_down, M, I_s)
    C = 1e-7 * 1e9;
    b = 0.5 * abs(x_max - x_min);
    l = 0.5 * abs(z_down - z_up);
    x0 = 0.5 * (x_min + x_max);
    z0 = 0.5 * (z_up + z_down);
    x1 = x - x0 + b;
    x2 = x - x0 - b;
    z1 = z0 - z - l;   % h 
    z2 = z0 - z + l;   % h + 2l
    
    %%% analytical solution 1 --- Equation (3-1-89) refered to Zhining Guan(2005) 
    A = 0.5 * log(((z1*z1 + x2*x2)*(z2*z2 + x1*x1))/((z1*z1 + x1*x1)*(z2*z2 + x2*x2)));
    F = (atan((x1)/(z1)) - atan((x2)/(z1))) - (atan((x1)/(z2)) - atan((x2)/(z2)));
    Hax = C * 2 * M * (sin(pi* I_s / 180) * A - cos(pi* I_s / 180) * F);
    Za = C * 2 * M * (cos(pi * I_s / 180) * A + sin(pi * I_s / 180) * F);
    delta_T = C * 2 * M * (cos(pi * (90 - 2 * I_s) / 180) * A - sin(pi * (90 - 2 * I_s) / 180) * F);
    %%% analytical solution 2 --- Equation (2.1 - 2.6) refered to Suang liu(2020) 
%     E = log(((z2*z2 + x1*x1)* (z1*z1 + x2*x2))/(((z1*z1 + x1*x1))*(z2*z2 + x2*x2)));
%     F1 = atan((2*b*z1)/(z1*z1 + (x1 - b)*(x1 - b) - b*b)) - atan((2*b*z2)/(z2*z2 +(x1 - b)*(x1 - b) - b*b));
%     Hax = C * 2 * M * (0.5 * sin(pi * I_s / 180) * E - cos (pi * I_s / 180)* F1);
%     Za = C * 2 * M * (0.5 * cos(pi * I_s / 180) * E + sin (pi * I_s / 180)* F1);
%     delta_T = Hax * cos(pi * I_s / 180) + Za * sin(pi * I_s / 180);  %sqrt(Hax*Hax + Za*Za);
end
