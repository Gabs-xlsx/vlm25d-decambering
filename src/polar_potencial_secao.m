function P2D = polar_potencial_secao(Malha, j, Nx)
% POLAR_POTENCIAL_SECAO  Polar potencial 2D da propria linha media discretizada.

    ini = (j-1)*Nx + 1;  fim = j*Nx;
    P = Malha.Paineis(ini:fim);

    c = zeros(Nx,1);
    for i = 1:Nx; c(i) = P(i).Corda; end
    X = [0; cumsum(c)] / sum(c);
    xv  = X(1:end-1) + 0.25*diff(X);
    xcp = X(1:end-1) + 0.75*diff(X);

    dz = zeros(Nx,1);
    for i = 1:Nx; dz(i) = -P(i).n(1) / P(i).n(3); end
    if isfield(Malha,'deflexao') && isstruct(Malha.deflexao)
        for i = 1:Nx
            if isfield(P(i),'is_control') && P(i).is_control && isfield(P(i),'nome_controle')
                nm = matlab.lang.makeValidName(P(i).nome_controle);
                if isfield(Malha.deflexao, nm)
                    d = deg2rad(Malha.deflexao.(nm));
                    if isfield(Malha,'deflexao_tipo') && strcmpi(Malha.deflexao_tipo,'assimetrico') && P(i).CP(2) < 0
                        d = -d;
                    end
                    dz(i) = tan(atan(dz(i)) + d) - tan(d);
                end
            end
        end
    end

    A = zeros(Nx,Nx);
    for a = 1:Nx
        for b = 1:Nx
            A(a,b) = 1 / (2*pi*(xcp(a) - xv(b)));
        end
    end

    resolve = @(al) A \ (al - dz);
    cl  = @(g) 2*sum(g);
    cm  = @(g) -2*sum(g .* (xv - 0.25));

    g0 = resolve(0);      g1 = resolve(0.01);
    cl0 = cl(g0);         cl1 = cl(g1);

    P2D.a0       = (cl1 - cl0) / 0.01;
    P2D.alpha_L0 = -cl0 / P2D.a0;
    P2D.Cm_c4    = cm(resolve(0.10));
    P2D.Cl       = @(a_rad) P2D.a0 * (a_rad - P2D.alpha_L0);

    P2D.alpha_L0_visc = P2D.alpha_L0;
    P2D.Cm_c4_visc    = P2D.Cm_c4;
end
