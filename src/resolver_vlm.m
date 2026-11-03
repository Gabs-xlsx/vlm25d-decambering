function [Resultados] = resolver_vlm(Malha, A, CondicaoVoo, Aeronave, varargin)
% RESOLVER_VLM  Resolve o sistema linear do VLM e integra as cargas.

    persistent avisou_arg5
    if ~isempty(varargin) && isempty(avisou_arg5)
        warning('resolver_vlm:argObsoleto', ...
            ['O 5o argumento (usar_inducao_vp) foi removido: agora existe um ' ...
             'unico modelo de forca (KJ + Trefftz). Sera ignorado.']);
        avisou_arg5 = true;
    end

    Paineis = Malha.Paineis;
    N   = numel(Paineis);
    rho = obter_campo(CondicaoVoo, 'rho', 1.225);
    V   = CondicaoVoo.V_inf;

    alpha = deg2rad(CondicaoVoo.alpha);
    beta  = deg2rad(obter_campo(CondicaoVoo, 'beta', 0));
    p_ad  = obter_campo(CondicaoVoo, 'p', 0);
    q_ad  = obter_campo(CondicaoVoo, 'q', 0);
    r_ad  = obter_campo(CondicaoVoo, 'r', 0);

    Omega = [p_ad*(2*V)/Aeronave.Bref, q_ad*(2*V)/Aeronave.Cref, r_ad*(2*V)/Aeronave.Bref];
    CG    = [Aeronave.Xcg, 0, Aeronave.Zcg];

    V_fs = V * [cos(alpha)*cos(beta), -sin(beta), sin(alpha)*cos(beta)];

    RHS = zeros(N, 1);
    for i = 1:N
        V_loc  = V_fs + cross(Omega, Paineis(i).CP - CG);
        RHS(i) = -dot(V_loc, Paineis(i).n);
    end

    Gamma = A \ RHS;

    q_din = 0.5 * rho * V^2;
    F_tot = [0 0 0];  M_tot = [0 0 0];
    for i = 1:N
        A_no = Paineis(i).P1 + 0.25*(Paineis(i).P2 - Paineis(i).P1);
        B_no = Paineis(i).P4 + 0.25*(Paineis(i).P3 - Paineis(i).P4);
        dl   = B_no - A_no;

        r_vp  = (A_no + B_no)/2 - CG;
        V_loc = V_fs + cross(Omega, r_vp);
        F_i   = rho * Gamma(i) * cross(V_loc, dl);

        Paineis(i).Gamma = Gamma(i);
        Paineis(i).Forca = F_i;
        Paineis(i).dCp   = 2*Gamma(i) / (V * Paineis(i).Corda);

        F_tot = F_tot + F_i;
        M_tot = M_tot + cross(r_vp, F_i);
    end

    Lift = F_tot(3)*cos(alpha) - F_tot(1)*sin(alpha);

    [w_ll, CDi, Gamma_faixa] = downwash_trefftz(Malha, Gamma, Aeronave, rho, V);

    Resultados.CL          = Lift / (q_din * Aeronave.Sref);
    Resultados.CDi         = CDi;
    Resultados.Cm          = M_tot(2) / (q_din * Aeronave.Sref * Aeronave.Cref);
    Resultados.Gamma       = Gamma;
    Resultados.Gamma_faixa = Gamma_faixa;
    Resultados.w_induzido  = w_ll;
    Resultados.alpha_ind   = rad2deg(atan2(w_ll, V));
    Resultados.Paineis     = Paineis;
end

function v = obter_campo(s, nome, padrao)
    if isfield(s, nome) && ~isempty(s.(nome)); v = s.(nome); else; v = padrao; end
end
