function [Cl, Cm, c_loc, Area] = cargas_por_faixa(Malha, Gamma, V_fs, rho, V)
% CARGAS_POR_FAIXA  Cl e Cm (em torno de c/4 local) de cada faixa de envergadura.

    Nx = Malha.Nx;  Ny = Malha.Ny;
    q  = 0.5*rho*V^2;
    Cl = zeros(Ny,1); Cm = zeros(Ny,1); c_loc = zeros(Ny,1); Area = zeros(Ny,1);

    for j = 1:Ny
        ini = (j-1)*Nx + 1;  fim = j*Nx;
        P = Malha.Paineis;

        P_ba = (P(ini).P1 + P(ini).P4)/2;
        P_bf = (P(fim).P2 + P(fim).P3)/2;
        x_c4 = P_ba(1) + 0.25*(P_bf(1) - P_ba(1));

        L = 0; M = 0; A = 0;
        for p = ini:fim
            A_no = P(p).P1 + 0.25*(P(p).P2 - P(p).P1);
            B_no = P(p).P4 + 0.25*(P(p).P3 - P(p).P4);
            Fz = rho * Gamma(p) * subsref(cross(V_fs, B_no - A_no), substruct('()',{3}));
            L  = L + Fz;
            M  = M - Fz * ((A_no(1)+B_no(1))/2 - x_c4);
            A  = A + P(p).Area;
        end
        env = norm(P(ini).P4 - P(ini).P1);
        c_loc(j) = A / env;   Area(j) = A;
        Cl(j) = L / (q*A);
        Cm(j) = M / (q*A*c_loc(j));
    end
end
