function [dCl, dCm] = derivadas_polares(Polares, Aeronave, Malha, Nx, a_ef_deg)
% DERIVADAS_POLARES  dCl_visc/dalpha_ef e dCm_visc/dalpha_ef [1/rad], por faixa.

    h = 2.0;
    [clp, cmp] = interpolar_polares(Polares, Aeronave, Malha, Nx, a_ef_deg + h);
    [clm, cmm] = interpolar_polares(Polares, Aeronave, Malha, Nx, a_ef_deg - h);
    dCl = (clp - clm) / (2*deg2rad(h));
    dCm = (cmp - cmm) / (2*deg2rad(h));
end
