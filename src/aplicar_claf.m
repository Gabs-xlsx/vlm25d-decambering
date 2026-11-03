function Malha = aplicar_claf(Malha, claf)
%APLICAR_CLAF  Corrige a inclinacao de sustentacao para a ESPESSURA do perfil.
%   Malha = aplicar_claf(Malha, 1.097)

    if nargin < 2 || isempty(claf) || abs(claf - 1) < 1e-12
        return;
    end
    if claf <= 0 || claf > 2
        error('aplicar_claf:valor', ...
              'claf = %g fora do plausivel (esperado ~1,00 a 1,20).', claf);
    end
    for i = 1:numel(Malha.Paineis)
        p = Malha.Paineis(i);
        Malha.Paineis(i).CP = p.VP + claf * (p.CP - p.VP);
    end
    Malha.claf = claf;
end
