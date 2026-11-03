function [Res, Malha_Def, hist] = resolver_vlm_25D(Malha, CasoVoo, Aeronave, Polares, opts)
% RESOLVER_VLM_25D  Acoplamento viscoso-invicido por decambering.

    if nargin < 5; opts = struct(); end
    tol      = obter(opts,'tol',1e-3);
    max_iter = obter(opts,'max_iter',200);
    Damp     = obter(opts,'damp',1.0);
    x2       = obter(opts,'x2',0.8);
    usar_tara= obter(opts,'usar_tara',true);
    verbose  = obter(opts,'verbose',true);

    modo_lib = obter(opts,'modo_liberacao','segmento');
    if strcmpi(modo_lib,'adaptativo') && ~obter(opts,'passe2',false)
        o1 = opts;
        o1.n_extrap = 0;  o1.n_borda = 0;  o1.livres = [];  o1.eta_lib = 0;
        o1.modo_liberacao = 'copia';  o1.passe2 = true;  o1.verbose = false;
        o1.max_iter = min(max_iter, obter(opts,'iter_sonda',25));
        R1 = resolver_vlm_25D(Malha, CasoVoo, Aeronave, Polares, o1);

        d1p = R1.delta1(:);  d2p = R1.delta2(:);
        djs = deflexao_por_faixa(Malha,'');  djs = djs(:);
        Nys = numel(d1p);
        segp = ones(Nys,1);
        for jj = 2:Nys; segp(jj) = segp(jj-1) + double(abs(djs(jj)-djs(jj-1)) > 1e-6); end
        osc = zeros(Nys,1);
        for jj = 1:Nys
            viz = [];
            if jj > 1   && segp(jj-1) == segp(jj); viz(end+1) = jj-1; end
            if jj < Nys && segp(jj+1) == segp(jj); viz(end+1) = jj+1; end
            if isempty(viz); continue; end
            osc(jj) = max(abs(d1p(jj) - mean(d1p(viz))), abs(d2p(jj) - mean(d2p(viz))));
        end
        nc = obter(opts,'cand_borda',4);
        cand = false(Nys,1);
        cand(1:min(nc,Nys)) = true;  cand(max(1,Nys-nc+1):Nys) = true;
        for jj = 1:Nys-1
            if abs(djs(jj) - djs(jj+1)) > 1e-6
                cand(max(1,jj-nc+1) : min(Nys,jj+nc)) = true;
            end
        end
        if any(~cand); base = median(osc(~cand)); else; base = 0; end
        if ~isfinite(base) || base < 0; base = 0; end
        lim_o = max(deg2rad(obter(opts,'osc_min',1)), 5*base);
        duro  = R1.clip1(:) | R1.clip2(:) | R1.fora_tabela(:);
        ruim  = duro | (cand & osc > lim_o);
        Nys   = numel(ruim);
        if sum(ruim) > 0.5*Nys; ruim = duro; end
        if sum(ruim) > 0.5*Nys
            error('vlm25D:sonda', ...
                  ['a sonda adaptativa reprovou %d de %d faixas. Isso nao e'' ' ...
                   'caso de liberar borda: a polar ou o caso de voo e'' que nao ' ...
                   'fecham.'], sum(ruim), Nys);
        end
        ruim = fechar_buracos_seg(ruim, deflexao_por_faixa(Malha,''));

        o2 = opts;
        o2.n_extrap = 0;  o2.n_borda = 0;  o2.livres = find(ruim);  o2.eta_lib = 0;
        o2.modo_liberacao = 'segmento';  o2.passe2 = true;
        [Res, Malha_Def, hist] = resolver_vlm_25D(Malha, CasoVoo, Aeronave, Polares, o2);
        Res.sonda = struct('n_liberadas',sum(ruim), 'residuo_sonda',R1.residuo, ...
                           'clip',sum(R1.clip1(:)|R1.clip2(:)), ...
                           'fora',sum(R1.fora_tabela(:)), 'faixas',find(ruim(:)).', ...
                           'osc_graus',rad2deg(osc(:)).', 'lim_osc_graus',rad2deg(lim_o), ...
                           'n_candidatas',sum(cand));
        return
    end
    lim_estado= obter(opts,'lim_estado',deg2rad(15));
    lim_estado= obter(opts,'lim_delta',lim_estado);
    lim_passo = obter(opts,'lim_passo', deg2rad(25));
    n_ponta  = obter(opts,'n_ponta',0);
    n_extrap = obter(opts,'n_extrap',0);
    esquema  = obter(opts,'esquema',2);
    estol_ant= obter(opts,'estol_inicial',[]);

    Nx = Malha.Nx;  Ny = Malha.Ny;
    eta_lib = obter(opts,'eta_lib',0);
    if eta_lib > 0 && n_extrap >= 0
        eta_j = abs(arrayfun(@(j) Malha.Paineis((j-1)*Nx+1).CP(2), 1:Ny)) / (Aeronave.Bref/2);
        dentro = find(eta_j <= eta_lib);
        if ~isempty(dentro)
            n_extrap = max([n_extrap, dentro(1)-1, Ny-dentro(end)]);
        end
    end
    dae1_d = nan(Ny,1);  dae2_d = nan(Ny,1);  dClv_d = nan(Ny,1);  dCmv_d = nan(Ny,1);
    prop_d = NaN;  J_cond = NaN;  J_svmin = NaN;  r_faixa = nan(Ny,1);
    a_t_v  = nan(Ny,1);  m_traj = nan(Ny,1);  n_raiz_v = zeros(Ny,1);
    modo_raiz = obter(opts,'modo_raiz','max');
    jac_cons = obter(opts,'jac_consistente', esquema == 2);
    rho = obter(CasoVoo,'rho',1.225);  V = CasoVoo.V_inf;
    al = deg2rad(CasoVoo.alpha);  be = deg2rad(obter(CasoVoo,'beta',0));
    V_fs = V * [cos(al)*cos(be), -sin(be), sin(al)*cos(be)];

    claf = obter(opts,'claf',1);
    Malha = aplicar_claf(Malha, claf);

    C = montar_tensor_influencia(Malha);
    D = decambering_dobradicas(Malha, x2);

    A0 = aic_de_normais(C, Malha.Paineis);
    [E_Cl, E_Cm, P2D] = tara_reducao_25d(Malha, A0, CasoVoo, Aeronave, D, Polares);
    if ~usar_tara; E_Cl(:) = 0; E_Cm(:) = 0; end
    if verbose
        fprintf('  Tara da reducao 3D->2D: |E_Cl| max = %.2e | |E_Cm| max = %.4f\n', ...
                max(abs(E_Cl)), max(abs(E_Cm)));
    end

    d1 = zeros(Ny,1);  d2 = zeros(Ny,1);
    if isempty(estol_ant); estol_ant = false(Ny,1); end

    ativas = true(Ny,1);
    ativas(1:n_ponta) = false;  ativas(end-n_ponta+1:end) = false;
    resolvidas = true(Ny,1);
    if n_extrap > 0
        resolvidas(1:n_extrap) = false;  resolvidas(end-n_extrap+1:end) = false;
    end
    livres = obter(opts,'livres',[]);
    if ~isempty(livres)
        resolvidas(livres) = false;
    end
    n_borda = obter(opts,'n_borda',1);
    if n_borda > 0
        dj_b = deflexao_por_faixa(Malha, '');
        for j = 1:Ny-1
            if abs(dj_b(j) - dj_b(j+1)) > 1e-6
                resolvidas(max(1,j-n_borda+1) : min(Ny,j+n_borda)) = false;
            end
        end
    end
    if ~any(resolvidas)
        error('vlm25D:sem_incognita', ...
              'todas as %d faixas foram liberadas -- nao sobrou incognita.', Ny);
    end
    idx_res = find(resolvidas);   M = numel(idx_res);
    yf_lib = zeros(Ny,1);
    for j = 1:Ny
        p1 = Malha.Paineis((j-1)*Nx + 1);
        yf_lib(j) = 0.5*(p1.P1(2) + p1.P4(2));
    end
    [S, D_lib] = mapa_liberacao(resolvidas, deflexao_por_faixa(Malha,''), ...
                                yf_lib, modo_lib, obter(opts,'t_extrap',2), ...
                                obter(opts,'n_ajuste',2));
    T = blkdiag(S, S);
    lin_sel = [idx_res(:); Ny + idx_res(:)];
    Malha_Def = Malha;
    hist = struct('iter',[],'erro',[],'CL',[],'Cm',[]);

    for it = 1:max_iter
        [Malha_Def.Paineis, dn1, dn2] = aplicar_decambering(Malha.Paineis, D, d1, d2, true);

        A   = aic_de_normais(C, Malha_Def.Paineis);
        n   = reshape([Malha_Def.Paineis.n], 3, Nx*Ny)';
        A_lu = decomposition(A, 'lu');
        Gam = A_lu \ (-n * V_fs(:));

        [Cl3, Cm3] = cargas_por_faixa(Malha_Def, Gam, V_fs, rho, V);

        a_ef = zeros(Ny,1);
        for j = 1:Ny
            a_ef(j) = Cl3(j)/P2D(j).a0 + P2D(j).alpha_L0 ...
                      - d1(j) - d2(j)*(1 - D.theta2/pi + sin(D.theta2)/pi);
        end
        a_ef_deg = rad2deg(a_ef);

        if esquema == 1
            [Cl_v, Cm_v] = interpolar_polares(Polares, Aeronave, Malha_Def, Nx, a_ef_deg, P2D);
        else
            pert = deg2rad(0.25);
            [Pp, ~, ~] = aplicar_decambering(Malha.Paineis, D, d1+pert, d2, false);
            Mp = Malha_Def;  Mp.Paineis = Pp;
            np = reshape([Pp.n], 3, Nx*Ny)';
            Gp = aic_de_normais(C, Pp) \ (-np * V_fs(:));
            Cl_p = cargas_por_faixa(Mp, Gp, V_fs, rho, V);
            a_p  = zeros(Ny,1);
            for j = 1:Ny
                a_p(j) = rad2deg(Cl_p(j)/P2D(j).a0 + P2D(j).alpha_L0 ...
                         - (d1(j)+pert) - d2(j)*(1 - D.theta2/pi + sin(D.theta2)/pi));
            end
            Cl_v = zeros(Ny,1);  Cm_v = zeros(Ny,1);
            for j = 1:Ny
                [PA, PCL, PCM] = polar_da_faixa(Polares, Aeronave, Malha_Def, Nx, j);
                [Cl_v(j), Cm_v(j), a_t_v(j), estol_ant(j), n_raiz_v(j)] = ...
                    alvo_esquema2(PA, PCL, PCM, a_ef_deg(j), Cl3(j), ...
                                  a_p(j), Cl_p(j), estol_ant(j), modo_raiz);
                den = a_p(j) - a_ef_deg(j);
                if abs(den) < 1e-9; den = sign(den + eps)*1e-9; end
                m_traj(j) = (Cl_p(j) - Cl3(j)) / den;
            end
        end
        R_Cl = Cl3 - Cl_v - E_Cl;
        R_Cm = Cm3 - Cm_v - E_Cm;

        erro = max(max(abs(R_Cl(resolvidas))), max(abs(R_Cm(resolvidas))));
        hist.iter(end+1) = it;  hist.erro(end+1) = erro;

        if verbose
            [eCl,jCl] = max(abs(R_Cl));  [eCm,jCm] = max(abs(R_Cm));
            fprintf(['  it %3d | erro %.5f  (Cl %.5f na faixa %2d | Cm %.5f na faixa %2d)' ...
                     ' | d1 raiz %+6.2f | d2 raiz %+7.2f | a_ef raiz %+6.2f deg\n'], ...
                    it, erro, eCl, jCl, eCm, jCm, rad2deg(d1(round(Ny/2))), ...
                    rad2deg(d2(round(Ny/2))), a_ef_deg(round(Ny/2)));
        end
        if erro < tol; break; end

        [dCl, dCm] = sensibilidades_decambering(Malha_Def, C, D, Gam, A_lu, V_fs, rho, V, dn1, dn2);

        if jac_cons && esquema == 2
            [dClv, dCmv] = derivadas_polares(Polares, Aeronave, Malha_Def, Nx, a_t_v);
            mr  = m_traj(:) * 180/pi;
            den_t = mr - dClv(:);
            den_t(abs(den_t) < 1e-6) = sign(den_t(abs(den_t) < 1e-6) + eps)*1e-6;
        else
            [dClv, dCmv] = derivadas_polares(Polares, Aeronave, Malha_Def, Nx, a_ef_deg);
            mr = nan(Ny,1);  den_t = nan(Ny,1);
        end
        dClv_d = dClv(:);  dCmv_d = dCmv(:);
        c2 = 1 - D.theta2/pi + sin(D.theta2)/pi;
        Jl1 = dCl{1};  Jl2 = dCl{2};  Jm1 = dCm{1};  Jm2 = dCm{2};
        for j = 1:Ny
            for k = 1:Ny
                dae1 = dCl{1}(j,k)/P2D(j).a0 - (j==k);
                dae2 = dCl{2}(j,k)/P2D(j).a0 - (j==k)*c2;
                if j == k; dae1_d(j) = dae1;  dae2_d(j) = dae2; end
                if jac_cons && esquema == 2
                    dat1 = (mr(j)*dae1 - dCl{1}(j,k)) / den_t(j);
                    dat2 = (mr(j)*dae2 - dCl{2}(j,k)) / den_t(j);
                    Jl1(j,k) = Jl1(j,k) - dClv(j)*dat1;
                    Jl2(j,k) = Jl2(j,k) - dClv(j)*dat2;
                    Jm1(j,k) = Jm1(j,k) - dCmv(j)*dat1;
                    Jm2(j,k) = Jm2(j,k) - dCmv(j)*dat2;
                else
                    Jl1(j,k) = Jl1(j,k) - dClv(j)*dae1;
                    Jl2(j,k) = Jl2(j,k) - dClv(j)*dae2;
                    Jm1(j,k) = Jm1(j,k) - dCmv(j)*dae1;
                    Jm2(j,k) = Jm2(j,k) - dCmv(j)*dae2;
                end
            end
        end
        J = [Jl1, Jl2; Jm1, Jm2];

        prop_d  = norm(dCl{2} - c2*dCl{1},'fro') / max(norm(dCl{2},'fro'),eps);
        sv      = svd(J);
        J_cond  = sv(1)/max(sv(end),eps);
        J_svmin = sv(end);
        r_faixa = dClv(:) ./ [P2D.a0]';

        idx = [~ativas; ~ativas];
        J(idx,:) = 0;  J(:,idx) = 0;
        J(sub2ind(size(J), find(idx), find(idx))) = 1;

        if n_extrap > 0
            J_red  = J(lin_sel, :) * T;
            r_red  = [R_Cl(idx_res); R_Cm(idx_res)];
            p_red  = -J_red \ r_red;
            p_red  = max(min(p_red, lim_passo), -lim_passo);
            passo  = T * p_red;
        else
            passo = -J \ [R_Cl; R_Cm];
            passo = max(min(passo, lim_passo), -lim_passo);
        end

        d1 = d1 + Damp*passo(1:Ny);
        d2 = d2 + Damp*passo(Ny+1:end);
        d1 = max(min(d1, lim_estado), -lim_estado);
        d2 = max(min(d2, lim_estado), -lim_estado);
    end

    Res = resolver_vlm(Malha_Def, aic_de_normais(C, Malha_Def.Paineis), CasoVoo, Aeronave);
    Res.delta1 = d1;  Res.delta2 = d2;  Res.alpha_efetivo = a_ef_deg;
    Res.E_Cl = E_Cl;  Res.E_Cm = E_Cm;  Res.residuo = erro;
    Res.convergiu = erro < tol;
    Res.estol = estol_ant;
    Res.esquema = esquema;
    Res.faixas_ativas = ativas;
    Res.n_ponta = n_ponta;
    Res.Cl_faixa = Cl3;  Res.R_Cl = R_Cl;  Res.R_Cm = R_Cm;
    Res.Cm_faixa = Cm3;
    Res.lim_estado = lim_estado;   Res.lim_passo = lim_passo;
    Res.lim_delta  = lim_estado;
    Res.clip1 = abs(d1) >= 0.999*lim_estado;
    Res.clip2 = abs(d2) >= 0.999*lim_estado;
    Res.c2 = 1 - D.theta2/pi + sin(D.theta2)/pi;
    Res.n_extrap = n_extrap;  Res.resolvidas = resolvidas;
    Res.lib = D_lib;
    Res.residuo_total = max(max(abs(R_Cl)), max(abs(R_Cm)));
    Res.alpha_alvo = a_t_v;   Res.m_traj = m_traj;
    Res.n_raiz = n_raiz_v;   Res.modo_raiz = modo_raiz;
    Res.max_raizes = max(n_raiz_v);
    Res.jac_consistente = jac_cons;
    Res.desloc_alvo = a_t_v - a_ef_deg;
    Res.prop_dCl = prop_d;
    Res.J_cond   = J_cond;   Res.J_svmin = J_svmin;
    Res.r_faixa  = r_faixa;
    Res.dae1 = dae1_d;  Res.dae2 = dae2_d;
    Res.dClv = dClv_d;  Res.dCmv = dCmv_d;
    if isfield(Malha_Def,'Pol_faixa') && ~isempty(Malha_Def.Pol_faixa) && ...
       isfield(Malha_Def.Pol_faixa,'alpha_min_real')
        amin_j = [Malha_Def.Pol_faixa.alpha_min_real]';
        amax_j = [Malha_Def.Pol_faixa.alpha_max_real]';
        Res.alpha_min_real = amin_j;  Res.alpha_max_real = amax_j;
        Res.fora_tabela = a_ef_deg(:) < amin_j | a_ef_deg(:) > amax_j;
        Res.n_fora = sum(Res.fora_tabela);
        if Res.n_fora > 0 && verbose
            fprintf(['  %d de %d faixas fora do dado real da SUA propria polar ' ...
                     '(alpha = %+.1f)\n'], Res.n_fora, Ny, CasoVoo.alpha);
        end
    elseif isfield(Polares,'alpha_min_real')
        amin = max([Polares.alpha_min_real]);  amax = min([Polares.alpha_max_real]);
        Res.alpha_min_real = amin;  Res.alpha_max_real = amax;
        Res.fora_tabela = a_ef_deg < amin | a_ef_deg > amax;
        Res.n_fora = sum(Res.fora_tabela);
        if Res.n_fora > 0
            warning('vlm25D:fora_tabela', ...
              ['%d de %d faixas operam FORA do dado real da polar [%+.1f %+.1f] ' ...
               'em alpha = %+.1f (alpha_ef de %+.1f a %+.1f). O alvo delas e ' ...
               'extrapolacao Viterna, nao medida.'], Res.n_fora, Ny, amin, amax, ...
              CasoVoo.alpha, min(a_ef_deg), max(a_ef_deg));
        end
    else
        Res.fora_tabela = false(Ny,1);  Res.n_fora = 0;
        Res.alpha_min_real = NaN;  Res.alpha_max_real = NaN;
    end
    if any(Res.clip1 | Res.clip2)
        warning('vlm25D:clip', ...
          ['delta saturou o limite FISICO de %d deg em alpha = %+.1f: ' ...
           '%d faixas em delta1, %d em delta2. Isso e resultado, nao erro -- ' ...
           'o perfil pede mais decambering do que uma dobradica representa.'], ...
          round(rad2deg(lim_estado)), CasoVoo.alpha, sum(Res.clip1), sum(Res.clip2));
    end
    if verbose
        if Res.convergiu
            fprintf('  >>> convergido em %d iteracoes (erro %.2e)\n', it, erro);
        else
            fprintf('  >>> NAO convergiu: %d iteracoes, erro %.4f\n', it, erro);
        end
    end
end

function v = obter(s,nome,padrao)
    if isfield(s,nome) && ~isempty(s.(nome)); v = s.(nome); else; v = padrao; end
end

function r = fechar_buracos_seg(r, dj)
    r = logical(r(:));  dj = dj(:);  n = numel(r);
    seg = ones(n,1);
    for j = 2:n; seg(j) = seg(j-1) + double(abs(dj(j)-dj(j-1)) > 1e-6); end
    for j = 2:n-1
        if r(j) || ~r(j-1) || ~r(j+1); continue; end
        if seg(j-1) ~= seg(j) || seg(j+1) ~= seg(j); continue; end
        if sum(~r(seg == seg(j))) <= 2; continue; end
        r(j) = true;
    end
end
