function [Cl, Cm] = interpolar_polares(Polares, Aeronave, Malha, Nx, a_ef_deg, P2D)
% INTERPOLAR_POLARES  Alvos viscosos por faixa, no alpha efetivo local.

    Ny = Malha.Ny;
    Cl = zeros(Ny,1);  Cm = zeros(Ny,1);

    if isfield(Malha,'Pol_faixa') && ~isempty(Malha.Pol_faixa)
        PF = Malha.Pol_faixa;
        for j = 1:Ny
            Cl(j) = interp1(PF(j).alpha, PF(j).Cl, a_ef_deg(j), 'linear','extrap');
            Cm(j) = interp1(PF(j).alpha, PF(j).Cm, a_ef_deg(j), 'linear','extrap');
        end
        return;
    end

    Sec = Aeronave.Asa.Secoes;

    for j = 1:Ny
        ini = (j-1)*Nx + 1;
        y = abs(Malha.Paineis(ini).CP(2));

        k = find(y >= [Sec(1:end-1).Y]-1e-6 & y <= [Sec(2:end).Y]+1e-6, 1, 'first');
        if isempty(k); k = numel(Sec)-1; end
        f = (y - Sec(k).Y) / (Sec(k+1).Y - Sec(k).Y);
        f = min(max(f,0),1);

        PR = Polares(min(k,numel(Polares)));
        PT = Polares(min(k+1,numel(Polares)));

        clR = interp1(PR.alpha, PR.Cl, a_ef_deg(j), 'linear', 'extrap');
        clT = interp1(PT.alpha, PT.Cl, a_ef_deg(j), 'linear', 'extrap');
        cmR = interp1(PR.alpha, PR.Cm, a_ef_deg(j), 'linear', 'extrap');
        cmT = interp1(PT.alpha, PT.Cm, a_ef_deg(j), 'linear', 'extrap');

        Cl(j) = clR + f*(clT - clR);
        Cm(j) = cmR + f*(cmT - cmR);

    end
end
