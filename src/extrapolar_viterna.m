function Polar_Ext = extrapolar_viterna(Polar, opts)
% EXTRAPOLAR_VITERNA  Estende a polar 2D para altos angulos (Viterna-Corrigan).

    if nargin < 2; opts = struct(); end
    B1        = obter(opts,'B1',1.11);
    xcp_fim   = obter(opts,'xcp_final',0.50);
    a_emenda  = obter(opts,'alpha_emenda',[]);

    a  = Polar.alpha(:);  Cl = Polar.Cl(:);
    Cd = Polar.Cd(:);     Cm = Polar.Cm(:);

    [~, i_max] = max(Cl);

    if isempty(a_emenda)
        i_e = escolher_emenda(a, Cl, Cd, i_max);
    else
        [~, i_e] = min(abs(a - a_emenda));
    end
    a_e = a(i_e);

    r_e   = deg2rad(a_e);
    Cl_e  = Cl(i_e);  Cd_e = Cd(i_e);
    B2 = (Cd_e - B1*sin(r_e)^2) / cos(r_e);
    A1 = B1/2;
    A2 = (Cl_e - Cd_e*sin(r_e)*cos(r_e)) * sin(r_e)/cos(r_e)^2;

    vit_Cl = @(r) A1*sin(2*r) + A2*(cos(r).^2 ./ max(sin(r),1e-6));
    vit_Cd = @(r) B1*sin(r).^2 + B2*cos(r);

    a_pos = (ceil(a_e)+1 : 1 : 90)';
    r_pos = deg2rad(a_pos);
    Cl_pos = vit_Cl(r_pos);
    Cd_pos = vit_Cd(r_pos);

    xcp_e   = xcp_do_dado(Cl_e, Cd_e, Cm(i_e), r_e);
    w       = min(max((a_pos - a_e)/(40 - a_e), 0), 1);
    xcp_pos = xcp_e + (xcp_fim - xcp_e).*w;
    Cm_pos  = -(Cl_pos.*cos(r_pos) + Cd_pos.*sin(r_pos)) .* (xcp_pos - 0.25);

    Cl_pos = suavizar(Cl_pos, a_pos, a_e, Cl_e - vit_Cl(deg2rad(a_e)));
    Cd_pos = suavizar(Cd_pos, a_pos, a_e, Cd_e - vit_Cd(deg2rad(a_e)));
    Cm_pos = suavizar(Cm_pos, a_pos, a_e, Cm(i_e) - Cm_pos(1));

    a_neg = []; Cl_neg = []; Cd_neg = []; Cm_neg = [];
    if a(1) > -30
        a_0  = a(1);  Cl_0 = Cl(1);  Cd_0 = Cd(1);  Cm_0 = Cm(1);
        m_cl = (Cl(2) - Cl(1)) / (a(2) - a(1));
        m_cm = (Cm(2) - Cm(1)) / (a(2) - a(1));

        a_neg = (-90 : 1 : floor(a_0)-1)';

        Cl_lin  = Cl_0 + m_cl*(a_neg - a_0);
        Cl_sat  = -0.65 * max(Cl);
        Cl_neg  = max(Cl_lin, Cl_sat);

        w_far   = min(max((-40 - a_neg)/50, 0), 1);
        Cl_neg  = Cl_neg .* (1 - w_far);

        Cd_neg  = Cd_0 + B1*sin(deg2rad(a_neg - a_0)).^2;

        Cm_neg  = Cm_0 + m_cm*(a_neg - a_0);
        Cm_neg  = Cm_neg .* (1 - w_far);

        Cl_neg = suavizar(Cl_neg, a_neg, a_0, Cl_0 - Cl_neg(end));
        Cd_neg = suavizar(Cd_neg, a_neg, a_0, Cd_0 - Cd_neg(end));
        Cm_neg = suavizar(Cm_neg, a_neg, a_0, Cm_0 - Cm_neg(end));
    end

    Polar_Ext        = Polar;
    Polar_Ext.alpha  = [a_neg;  a(1:i_e);  a_pos];
    Polar_Ext.Cl     = [Cl_neg; Cl(1:i_e); Cl_pos];
    Polar_Ext.Cd     = [Cd_neg; Cd(1:i_e); Cd_pos];
    Polar_Ext.Cm     = [Cm_neg; Cm(1:i_e); Cm_pos];
    Polar_Ext.alpha_emenda = a_e;
    Polar_Ext.alpha_stall  = a(i_max);
    Polar_Ext.alpha_min_real = a(1);
    Polar_Ext.alpha_max_real = a(i_e);
    Polar_Ext.Origem = sprintf('%s + Viterna (emenda em %.1f deg)', Polar.Origem, a_e);
end

function i_e = escolher_emenda(a, Cl, Cd, i_max)
    n = numel(a);
    i_e = i_max;
    for i = i_max+1 : n
        salto_cd = Cd(i) > 2.5*Cd(i-1);
        salto_cl = abs(Cl(i) - Cl(i-1)) > 0.25;
        if salto_cd || salto_cl; break; end
        i_e = i;
    end
end

function xcp = xcp_do_dado(Cl, Cd, Cm, r)
    N = Cl*cos(r) + Cd*sin(r);
    if abs(N) < 1e-6; xcp = 0.25; else; xcp = 0.25 - Cm/N; end
    xcp = min(max(xcp, 0.20), 0.60);
end

function y = suavizar(y, a, a_e, offset)
    if abs(offset) < 1e-12; return; end
    w = max(1 - abs(a - a_e)/10, 0);
    y = y + offset*w;
end

function v = obter(s,nome,padrao)
    if isfield(s,nome) && ~isempty(s.(nome)); v = s.(nome); else; v = padrao; end
end
