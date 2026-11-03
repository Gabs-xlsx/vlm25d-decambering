function [PA, PCL, PCM] = polar_da_faixa(Polares, Aeronave, Malha, Nx, j)
% POLAR_DA_FAIXA  Polar interpolada na envergadura para a faixa j, como CURVA

    if isfield(Malha,'Pol_faixa') && ~isempty(Malha.Pol_faixa)
        PF  = Malha.Pol_faixa(j);
        PA  = PF.alpha(:);  PCL = PF.Cl(:);  PCM = PF.Cm(:);
        return;
    end

    Sec = Aeronave.Asa.Secoes;
    y = abs(Malha.Paineis((j-1)*Nx+1).CP(2));

    k = find(y >= [Sec(1:end-1).Y]-1e-6 & y <= [Sec(2:end).Y]+1e-6, 1, 'first');
    if isempty(k); k = numel(Sec)-1; end
    f = min(max((y - Sec(k).Y)/(Sec(k+1).Y - Sec(k).Y), 0), 1);

    PR = Polares(min(k,numel(Polares)));
    PT = Polares(min(k+1,numel(Polares)));

    PA  = unique([PR.alpha(:); PT.alpha(:)]);
    clR = interp1(PR.alpha, PR.Cl, PA, 'linear', 'extrap');
    clT = interp1(PT.alpha, PT.Cl, PA, 'linear', 'extrap');
    cmR = interp1(PR.alpha, PR.Cm, PA, 'linear', 'extrap');
    cmT = interp1(PT.alpha, PT.Cm, PA, 'linear', 'extrap');

    PCL = clR + f*(clT - clR);
    PCM = cmR + f*(cmT - cmR);
end
