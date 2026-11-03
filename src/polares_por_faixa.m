function [Malha, Rel] = polares_por_faixa(Malha, Polares, Aeronave, Fam, opts)
%POLARES_POR_FAIXA  Monta a polar-alvo de CADA faixa: interpola na envergadura,
%   [Malha, Rel] = polares_por_faixa(Malha, Polares, Aeronave, Fam)

    if nargin < 4; Fam = []; end
    if nargin < 5; opts = struct(); end
    nome   = obter(opts,'nome','');
    modo   = lower(obter(opts,'se_faltar','erro'));
    tol    = obter(opts,'tol',0.25);
    verb   = obter(opts,'verbose',false);

    x_dob = obter(opts,'x_dob',[]);
    if isempty(x_dob)
        x_dob = 0.80;
        if isfield(Aeronave.Asa,'Controles') && ~isempty(Aeronave.Asa.Controles)
            x_dob = 1 - mean([Aeronave.Asa.Controles.corda_c]);
        end
    end

    Nx = Malha.Nx;  Ny = Malha.Ny;
    Sec = Aeronave.Asa.Secoes;

    [dj, info] = deflexao_por_faixa(Malha, nome);

    dfam = [];
    if ~isempty(Fam); dfam = [Fam.delta]; end

    PF = repmat(struct('alpha',[],'Cl',[],'Cm',[],'delta',0,'fam',0, ...
                       'y',0,'alpha_L0',NaN,'Cm0',NaN,'origem','', ...
                       'alpha_min_real',NaN,'alpha_max_real',NaN), 1, Ny);
    usou_teoria = false;  usou_limpo = false;

    for j = 1:Ny
        y = abs(Malha.Paineis((j-1)*Nx+1).CP(2));
        d = dj(j);

        if abs(d) < tol
            Pj = Polares;  fam = 0;  org = 'limpa';
        else
            m = [];
            if ~isempty(dfam)
                [e, m0] = min(abs(dfam - d));
                if e <= tol; m = m0; end
            end
            if ~isempty(m)
                Pj = Fam(m).Polares;  fam = m;
                org = obter(Fam(m),'origem','defletida');
            else
                switch modo
                    case 'teoria'
                        Pj = familia_teorica(Polares, d, x_dob);
                        fam = -1;  org = 'teoria';  usou_teoria = true;
                    case 'limpo'
                        Pj = Polares;  fam = 0;  org = 'LIMPA (faltou familia)';
                        usou_limpo = true;
                    otherwise
                        error('polares_por_faixa:sem_familia', ...
                          ['faixa %d (y/(b/2) = %+.3f) pede delta = %+.2f deg e ' ...
                           'nenhuma familia cobre esse valor (tol %.2f deg).\n' ...
                           'Familias disponiveis: %s\n' ...
                           'Gere a familia que falta ou passe opts.se_faltar = ''teoria''.'], ...
                          j, Malha.Paineis((j-1)*Nx+1).CP(2)/(Aeronave.Bref/2), ...
                          d, tol, mat2str(round(dfam,2)));
                end
            end
        end

        k = find(y >= [Sec(1:end-1).Y]-1e-6 & y <= [Sec(2:end).Y]+1e-6, 1, 'first');
        if isempty(k); k = numel(Sec)-1; end
        f = min(max((y - Sec(k).Y)/(Sec(k+1).Y - Sec(k).Y), 0), 1);

        PR = Pj(min(k,   numel(Pj)));
        PT = Pj(min(k+1, numel(Pj)));

        PA  = unique([PR.alpha(:); PT.alpha(:)]);
        clR = interp1(PR.alpha, PR.Cl, PA, 'linear','extrap');
        clT = interp1(PT.alpha, PT.Cl, PA, 'linear','extrap');
        cmR = interp1(PR.alpha, PR.Cm, PA, 'linear','extrap');
        cmT = interp1(PT.alpha, PT.Cm, PA, 'linear','extrap');

        PF(j).alpha = PA;
        PF(j).Cl    = clR + f*(clT - clR);
        PF(j).Cm    = cmR + f*(cmT - cmR);
        PF(j).delta = d;
        PF(j).fam   = fam;
        PF(j).origem= org;
        PF(j).y     = Malha.Paineis((j-1)*Nx+1).CP(2)/(Aeronave.Bref/2);
        PF(j).alpha_L0 = alpha_L0_de(PA, PF(j).Cl);
        PF(j).Cm0   = interp1(PA, PF(j).Cm, 0, 'linear','extrap');
        if isfield(PR,'alpha_min_real') && isfield(PT,'alpha_min_real') && ...
           ~isempty(PR.alpha_min_real) && ~isempty(PT.alpha_min_real)
            PF(j).alpha_min_real = max(PR.alpha_min_real, PT.alpha_min_real);
            PF(j).alpha_max_real = min(PR.alpha_max_real, PT.alpha_max_real);
        else
            PF(j).alpha_min_real = PA(1);   PF(j).alpha_max_real = PA(end);
        end
    end

    Malha.Pol_faixa = PF;
    Malha.por_faixa = true;

    Rel = struct();
    Rel.y        = [PF.y]';
    Rel.delta    = [PF.delta]';
    Rel.fam      = [PF.fam]';
    Rel.alpha_L0 = [PF.alpha_L0]';
    Rel.Cm0      = [PF.Cm0]';
    Rel.origem   = {PF.origem}';
    Rel.deltas   = unique(round(Rel.delta*1e3)/1e3);
    Rel.info_deflexao = info;
    Rel.x_dob    = x_dob;
    Rel.usou_teoria = usou_teoria;
    Rel.usou_limpo  = usou_limpo;

    if usou_limpo
        warning('polares_por_faixa:limpo', ...
          ['ALGUMA faixa defletida ficou com a polar LIMPA (se_faltar = ''limpo''). ' ...
           'O resultado NAO representa a asa com comando.']);
    end
    if verb
        fprintf('  seletor: %d faixas | deltas locais %s | dobradica x/c %.3f\n', ...
                Ny, mat2str(round(Rel.deltas',2)), x_dob);
        lim  = Rel.delta == 0;
        baixo = Rel.delta > 0;   cima = Rel.delta < 0;
        if any(lim)
            fprintf('  alpha_L0: %d limpas %+.2f deg', sum(lim), mean(Rel.alpha_L0(lim),'omitnan'));
            if any(baixo)
                fprintf(' | %d para baixo %+.2f (salto %+.2f)', sum(baixo), ...
                    mean(Rel.alpha_L0(baixo),'omitnan'), ...
                    mean(Rel.alpha_L0(baixo),'omitnan') - mean(Rel.alpha_L0(lim),'omitnan'));
            end
            if any(cima)
                fprintf(' | %d para cima %+.2f (salto %+.2f)', sum(cima), ...
                    mean(Rel.alpha_L0(cima),'omitnan'), ...
                    mean(Rel.alpha_L0(cima),'omitnan') - mean(Rel.alpha_L0(lim),'omitnan'));
            end
            fprintf('\n');
        end
    end
end

function a0 = alpha_L0_de(a, cl)
    a0 = NaN;
    s = sign(cl(:));
    k = find(s(1:end-1) < 0 & s(2:end) >= 0);
    if isempty(k); return; end
    [~,q] = min(abs(a(k) - 0));  i = k(q);
    if cl(i+1) == cl(i); a0 = a(i); return; end
    a0 = a(i) - cl(i)*(a(i+1)-a(i))/(cl(i+1)-cl(i));
end

function v = obter(s,nome,padrao)
    if isstruct(s) && isfield(s,nome) && ~isempty(s.(nome)); v = s.(nome); else; v = padrao; end
end
