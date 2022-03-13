%surf_integral
%第一类曲面积分
%   I = surf_integral(f, z, [x,y], [x_m,x_M], [y_m,y_M])
%   I = surf_integral(f, [x,y,z], [u,v], [u_m,u_M], [v_m,v_M])
%  Examples:
%  计算int_int(x^2*y+z*y^2)dS, 积分曲面如下：
%  x=ucosv, y=usinv, z=v, 0<=u<=a, 0<=v<=2*pi
%  MATLAB求解语句
%  syms u v; syms a positive;
%  x=u*cos(v); y=u*sin(v); z=v; f=x^2*y+z*y^2;
%  I = surf_integral(f,[x,y,z],[u,v],[0,a],[0,2*pi])
function I = surf_integral(f,vars,t,a,b)
if length(f)==1
    if length(vars)~=1
        E = simplify(sum(diff(vars,t(1)).^2));
        F = sum(diff(vars,t(1)).*diff(vars,t(2)));
        G = simplify(sum(diff(vars,t(2)).^2));
    else
        E = simplify(1+diff(vars,t(1))^2);
        F = diff(vars,t(1))*diff(vars,t(2));
        G = simplify(1+diff(vars,t(2))^2);
    end
    I = int(int(simplify(f*sqrt(E*G-F^2)),t(1),a(1),a(2)),t(2),b(1),b(2));
else
    if length(vars)~=1
        A = diff(vars(2),t(1))*diff(vars(3),t(2)) - diff(vars(3),t(1))*diff(vars(2),t(2));
        B = diff(vars(3),t(1))*diff(vars(1),t(2)) - diff(vars(1),t(1))*diff(vars(3),t(2));
        C = diff(vars(1),t(1))*diff(vars(2),t(2)) - diff(vars(2),t(1))*diff(vars(1),t(2));
    else
        A = - diff(vars,t(1));
        B = - diff(vars,t(2));
        C = 1;
    end
    f = f(:); abc = [A, B, C];
    I = int(int(simplify(abc*f),t(1),a(1),a(2)),t(2),b(1),b(2));
end
