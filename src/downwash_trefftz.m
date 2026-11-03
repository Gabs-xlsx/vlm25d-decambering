function [w_ll, CDi, Gamma_faixa] = downwash_trefftz(Malha, Gamma, Aeronave, rho, V)
% DOWNWASH_TREFFTZ  Angulo induzido e arrasto induzido pelo plano de Trefftz.

    Nx = Malha.Nx;  Ny = Malha.Ny;

    Gamma_faixa = sum(reshape(Gamma(:), Nx, Ny), 1)';
    y_no  = Malha.Y_nodes(:);
    y_c   = (y_no(1:end-1) + y_no(2:end)) / 2;
    dy    = diff(y_no);

    dGamma = diff([0; Gamma_faixa; 0]);

    w_ll = zeros(Ny, 1);
    for j = 1:Ny
        w_ll(j) = -(1/(4*pi)) * sum(dGamma ./ (y_c(j) - y_no));
    end

    Di  = -rho * sum(Gamma_faixa .* w_ll .* dy);
    CDi = Di / (0.5*rho*V^2 * Aeronave.Sref);
end
