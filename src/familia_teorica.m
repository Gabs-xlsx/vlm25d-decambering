function Pf = familia_teorica(Polares, delta_deg, x_dob)
%FAMILIA_TEORICA  Aplica polar_deslocada a TODAS as secoes de uma familia.
%   Pf = familia_teorica(Pol, 10, 0.80)

    n = numel(Polares);
    Pf = polar_deslocada(Polares(1), delta_deg, x_dob);
    if n > 1
        Pf = repmat(Pf, 1, n);
        for k = 2:n
            Pf(k) = polar_deslocada(Polares(k), delta_deg, x_dob);
        end
    end
end
