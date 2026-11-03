function ref = referencias_da_geometria(Aeronave)
% REFERENCIAS_DA_GEOMETRIA  Calcula Sref, Bref e Cref A PARTIR das secoes.

    S = Aeronave.Asa.Secoes;
    y = [S.Y];  c = [S.Corda];

    meia_asa = min(y) >= -1e-9;

    area = 0;  integral_c2 = 0;
    for k = 1:numel(y)-1
        dy = y(k+1) - y(k);
        area        = area + 0.5*(c(k) + c(k+1))*dy;
        integral_c2 = integral_c2 + (c(k)^2 + c(k)*c(k+1) + c(k+1)^2)*dy/3;
    end

    if meia_asa
        ref.Sref = 2*area;
        ref.Bref = 2*(max(y) - min(y));
        integral_c2 = 2*integral_c2;
    else
        ref.Sref = area;
        ref.Bref = max(y) - min(y);
    end
    ref.Cref = integral_c2 / ref.Sref;
    ref.AR   = ref.Bref^2 / ref.Sref;
    ref.meia_asa = meia_asa;
end
