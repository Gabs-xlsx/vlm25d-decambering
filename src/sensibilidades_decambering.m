function [dCl, dCm] = sensibilidades_decambering(Malha, C, D, Gamma, A_lu, V_fs, rho, V, dn_d1, dn_d2)
% SENSIBILIDADES_DECAMBERING  Jacobiana potencial 2N x 2N EXATA.

    Ny = numel(D.idx_faixa);
    dCl = {zeros(Ny,Ny), zeros(Ny,Ny)};
    dCm = {zeros(Ny,Ny), zeros(Ny,Ny)};

    for var = 1:2
        if var == 1; dn_all = dn_d1; else; dn_all = dn_d2; end
        for j = 1:Ny
            dn = dn_all(:,:,j);
            if ~any(dn(:)); continue; end

            dA = C(:,:,1).*dn(:,1) + C(:,:,2).*dn(:,2) + C(:,:,3).*dn(:,3);
            db = -dn * V_fs(:);

            dGamma = A_lu \ (db - dA*Gamma);

            [dcl, dcm] = cargas_por_faixa(Malha, dGamma, V_fs, rho, V);
            dCl{var}(:,j) = dcl;
            dCm{var}(:,j) = dcm;
        end
    end
end
