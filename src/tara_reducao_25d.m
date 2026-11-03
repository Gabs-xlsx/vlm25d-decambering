function [E_Cl, E_Cm, P2D] = tara_reducao_25d(Malha, A, CasoVoo, Aeronave, D, Polares)
% TARA_REDUCAO_25D  Vies da reducao 3D -> 2D, por faixa.

    Nx = Malha.Nx;  Ny = Malha.Ny;
    rho = obter(CasoVoo,'rho',1.225);  V = CasoVoo.V_inf;
    al  = deg2rad(CasoVoo.alpha);  be = deg2rad(obter(CasoVoo,'beta',0));
    V_fs = V * [cos(al)*cos(be), -sin(be), sin(al)*cos(be)];

    n = reshape([Malha.Paineis.n], 3, Nx*Ny)';
    Gamma = A \ (-n * V_fs(:));
    [Cl3, Cm3] = cargas_por_faixa(Malha, Gamma, V_fs, rho, V);

    E_Cl = zeros(Ny,1);  E_Cm = zeros(Ny,1);
    P2D  = polar_potencial_secao(Malha, 1, Nx);
    P2D  = repmat(P2D, Ny, 1);
    for j = 1:Ny
        P2D(j) = polar_potencial_secao(Malha, j, Nx);

        if nargin >= 6 && ~isempty(Polares)
            [PA, PCL, PCM] = polar_da_faixa(Polares, Aeronave, Malha, Nx, j);
            [~, iz] = min(abs(PCL));
            ac  = PA(iz);
            lin = PA >= ac - 4 & PA <= ac + 4;
            if sum(lin) < 3
                [~, ord] = sort(abs(PA - ac));
                lin = false(size(PA));
                lin(ord(1:min(5, numel(PA)))) = true;
            end
            cf  = polyfit(deg2rad(PA(lin)), PCL(lin), 1);
            P2D(j).alpha_L0_visc = -cf(2)/cf(1);
            P2D(j).Cm_c4_visc    = mean(PCM(lin));
        else
            P2D(j).alpha_L0_visc = P2D(j).alpha_L0;
            P2D(j).Cm_c4_visc    = P2D(j).Cm_c4;
        end
        a_ef   = Cl3(j)/P2D(j).a0 + P2D(j).alpha_L0;
        E_Cl(j) = Cl3(j) - P2D(j).Cl(a_ef);
        E_Cm(j) = Cm3(j) - P2D(j).Cm_c4;
    end
end

function v = obter(s,nome,padrao)
    if isfield(s,nome) && ~isempty(s.(nome)); v = s.(nome); else; v = padrao; end
end
