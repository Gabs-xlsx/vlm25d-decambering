function [Paineis, dn_d1, dn_d2] = aplicar_decambering(Paineis, D, d1, d2, quer_deriv)
% APLICAR_DECAMBERING  Gira as normais segundo delta1/delta2 e devolve as

    if nargin < 5; quer_deriv = false; end
    N = numel(Paineis);
    dn_d1 = [];  dn_d2 = [];
    if quer_deriv
        Ny = numel(D.idx_faixa);
        dn_d1 = zeros(N, 3, Ny);
        dn_d2 = zeros(N, 3, Ny);
    end

    for j = 1:numel(D.idx_faixa)
        for p = D.idx_faixa{j}
            th = d1(j);
            if any(D.idx_d2{j} == p); th = th + d2(j); end

            n0 = Paineis(p).n0;
            kx = cross(D.k_rot(p,:), n0);
            Paineis(p).n = n0*cos(th) + kx*sin(th);

            if quer_deriv
                dn = -n0*sin(th) + kx*cos(th);
                dn_d1(p,:,j) = dn;
                if any(D.idx_d2{j} == p); dn_d2(p,:,j) = dn; end
            end
        end
    end
end
