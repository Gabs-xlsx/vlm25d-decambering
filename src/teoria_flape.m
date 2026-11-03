function T = teoria_flape(x_dob)
%TEORIA_FLAPE  Coeficientes de flape por teoria de perfil fino, num so' lugar.
%   T = teoria_flape(0.80)

    if nargin < 1 || isempty(x_dob); x_dob = 0.80; end
    if x_dob <= 0 || x_dob >= 1
        error('teoria_flape:x','x_dob deve estar em (0,1); recebi %g.', x_dob);
    end

    th = acos(1 - 2*x_dob);

    T.x_dob = x_dob;
    T.theta = th;
    T.c2    = 1 - th/pi + sin(th)/pi;
    T.m2    = 0.25*sin(2*th) - 0.5*sin(th);
    T.dCl_ddelta = 2*pi*T.c2;
    T.dCm_ddelta = T.m2;

    T.dalphaL0_ddelta_deg = -T.c2;
    T.dCm_ddelta_deg      = T.m2 * pi/180;
    T.dCl_ddelta_deg      = T.dCl_ddelta * pi/180;

    T.dalphaL0 = @(d_deg) -T.c2 * rad2deg(tan(deg2rad(d_deg)));
end
