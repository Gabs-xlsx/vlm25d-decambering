function D = decambering_dobradicas(Malha, x2)
% DECAMBERING_DOBRADICAS  Geometria das duas funcoes de decambering por faixa.

    if nargin < 2 || isempty(x2); x2 = 0.8; end

    Nx = Malha.Nx;  Ny = Malha.Ny;

    D.x2     = x2;
    D.theta2 = acos(1 - 2*x2);

    D.dCl_d1 = 2*pi;
    D.dCm_d1 = 0;

    D.dCl_d2 = 2*(pi - D.theta2) + 2*sin(D.theta2);
    D.dCm_d2 = 0.25*sin(2*D.theta2) - 0.5*sin(D.theta2);

    D.idx_faixa = cell(Ny,1);
    D.idx_d2    = cell(Ny,1);
    for j = 1:Ny
        ini = (j-1)*Nx + 1;  fim = j*Nx;
        D.idx_faixa{j} = ini:fim;

        c_tot = sum([Malha.Paineis(ini:fim).Corda]);
        s = 0;  corte = fim + 1;
        for p = ini:fim
            s = s + Malha.Paineis(p).Corda;
            if s/c_tot >= x2; corte = p; break; end
        end
        if corte > fim; corte = fim; end
        D.idx_d2{j} = corte:fim;
    end

    D.k_rot = zeros(numel(Malha.Paineis), 3);
    for p = 1:numel(Malha.Paineis)
        k = cross(Malha.Paineis(p).n0, [1 0 0]);
        D.k_rot(p,:) = k / norm(k);
    end
end
