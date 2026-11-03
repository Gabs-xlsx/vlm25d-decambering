function Polar = run_xfoil(arquivo_perfil, Re, alpha_array, opts)
% RUN_XFOIL  Polar viscosa 2D pelo XFOIL, com varredura ancorada em DOIS processos.

    if nargin < 4; opts = struct(); end
    n_paineis = obter(opts,'n_paineis',240);
    dens_ba   = obter(opts,'dens_ba',1.5);
    raz_bf    = obter(opts,'raz_bf',0.15);
    if isfield(opts,'ancora'); anc = opts.ancora; else; anc = 4.0; end

    caminho_perfil = fullfile('Perfis', arquivo_perfil);
    if ~exist(caminho_perfil,'file')
        error('Perfil "%s" nao encontrado na pasta "Perfis".', arquivo_perfil);
    end

    aa = sort(alpha_array(:));
    cfg = struct('perfil',caminho_perfil,'Re',Re,'n',n_paineis, ...
                 'p',dens_ba,'t',raz_bf, ...
                 'viscoso', obter(opts,'viscoso',true));
    cfg.xtr   = obter(opts,'xtr',[1 1]);
    cfg.ncrit = obter(opts,'ncrit',9);
    cfg.limite_s = obter(opts,'limite_s',300);

    if isempty(anc)
        D = uma_varredura(cfg, aa, 1);
    else
        [~,ia] = min(abs(aa - anc));
        Dsobe  = uma_varredura(cfg, aa(ia:end),   2);
        Ddesce = uma_varredura(cfg, aa(ia:-1:1),  3);
        D = [Dsobe; Ddesce];
    end

    Polar = struct('Nome',arquivo_perfil,'Origem','XFOIL','Re',Re, ...
                   'alpha',[],'Cl',[],'Cd',[],'Cm',[], ...
                   'alpha_pedido',aa,'alpha_faltando',[]);

    if isempty(D)
        warning('xfoil:vazio','Polar do XFOIL vazia para %s (nada convergiu).', arquivo_perfil);
        Polar.alpha_faltando = aa;  return;
    end

    [~, iu] = unique(round(D(:,1)*100)/100, 'stable');
    D = D(iu,:);
    [~, io] = sort(D(:,1));
    D = D(io,:);

    Polar.alpha = D(:,1);  Polar.Cl = D(:,2);
    Polar.Cd    = D(:,3);  Polar.Cm = D(:,5);

    Polar.alpha_faltando = setdiff(round(aa*2)/2, round(Polar.alpha*2)/2);
    nf = numel(Polar.alpha_faltando);
    if nf > 0
        warning('xfoil:buracos', ...
          '%s Re %.2e: %d de %d alphas NAO convergiram (de %+.1f a %+.1f)', ...
          arquivo_perfil, Re, nf, numel(aa), min(Polar.alpha_faltando), ...
          max(Polar.alpha_faltando));
    end
end

function D = uma_varredura(cfg, alphas, tag)
    cmd = sprintf('xfoil_cmd_%d.txt', tag);
    pol = sprintf('xfoil_pol_%d.txt', tag);
    log = sprintf('xfoil_log_%d.txt', tag);
    for f = {cmd,pol,log}
        if exist(f{1},'file'); delete(f{1}); end
    end

    fid = fopen(cmd,'w');
    fprintf(fid,'LOAD %s\n', cfg.perfil);
    fprintf(fid,'PLOP\nG F\n\n');
    fprintf(fid,'PPAR\n');
    fprintf(fid,'N %d\n', cfg.n);
    fprintf(fid,'P %.2f\n', cfg.p);
    fprintf(fid,'T %.2f\n', cfg.t);
    fprintf(fid,'\n\n');
    fprintf(fid,'OPER\n');
    if ~isfield(cfg,'viscoso') || cfg.viscoso
        fprintf(fid,'VISC %f\n', cfg.Re);
        if any(abs(cfg.xtr(:).' - [1 1]) > 1e-9) || abs(cfg.ncrit - 9) > 1e-9
            fprintf(fid,'VPAR\nN %g\nXTR %.4f %.4f\n\n', cfg.ncrit, cfg.xtr(1), cfg.xtr(2));
        end
    end
    fprintf(fid,'ITER 200\n');
    fprintf(fid,'PACC\n');
    fprintf(fid,'%s\n\n', pol);
    for i = 1:numel(alphas)
        fprintf(fid,'ALFA %.2f\n', alphas(i));
    end
    fprintf(fid,'PACC\n');
    fprintf(fid,'\n\n\nQUIT\n\n');
    fclose(fid);

    [st, expirou] = xfoil_exec(cmd, log, obter(cfg,'limite_s',300));
    if expirou
        warning('xfoil:limite', ...
                'XFOIL passou de %g s e foi morto; os alphas restantes desta varredura viram buraco.', ...
                obter(cfg,'limite_s',300));
    elseif st ~= 0
        warning('xfoil:status','XFOIL retornou status %d (ver %s).', st, log);
    end

    D = [];
    if exist(pol,'file')
        f = fopen(pol,'r');
        for i = 1:12; fgetl(f); end
        D = fscanf(f,'%f %f %f %f %f %f %f',[7,inf])';
        fclose(f);
    end
    for f = {cmd,pol,log}
        if exist(f{1},'file'); delete(f{1}); end
    end
end

function v = obter(s,nome,padrao)
    if isfield(s,nome) && ~isempty(s.(nome)); v = s.(nome); else; v = padrao; end
end
