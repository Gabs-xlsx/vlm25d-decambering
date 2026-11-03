function [Pol, info] = polar_estoque(perfil, Re, alphas, opts)
% POLAR_ESTOQUE  Polar viscosa do XFOIL com estoque em disco e recuperacao
%   [Pol, info] = polar_estoque('s1223.dat', 2.44e5, -12:0.5:20);
%   [Pol, info] = polar_estoque('E423.dat', 2.44e5, -12:0.5:20, struct('forcar',true));

    if nargin < 4; opts = struct(); end
    pasta   = obter(opts,'pasta','Polares');
    forcar  = obter(opts,'forcar',false);
    recup   = obter(opts,'recuperar', obter(opts,'rampa',true));
    visc    = obter(opts,'viscoso',true);
    npan    = obter(opts,'n_paineis',240);
    anc     = obter(opts,'ancora',4);
    passo_r = obter(opts,'passo_rec',0.25);
    iter_r  = obter(opts,'iter_rec',60);
    verb    = obter(opts,'verbose',true);
    dba     = obter(opts,'dens_ba',1.5);
    rbf     = obter(opts,'raz_bf',0.15);
    xtr     = obter(opts,'xtr',[1 1]);
    ncr     = obter(opts,'ncrit',9);

    arq_dat = fullfile('Perfis', perfil);
    if exist(arq_dat,'file') ~= 2
        error('polar_estoque:perfil','"%s" nao encontrado em Perfis/.', perfil);
    end
    if exist(pasta,'dir') ~= 7; mkdir(pasta); end

    aa  = sort(alphas(:)).';
    Rer = str2double(sprintf('%.3g', Re));
    [~, base] = fileparts(perfil);
    hd  = md5_arquivo(arq_dat);
    ha  = md5_texto(sprintf('%.4f,', aa));
    modo = 'visc';  if ~visc; modo = 'inv'; end
    chave = sprintf('%s_%s_Re%.3g_a%s_N%d_%s', base, hd(1:8), Rer, ha(1:6), npan, modo);
    if any(abs(xtr(:).' - [1 1]) > 1e-9) || abs(ncr - 9) > 1e-9
        chave = sprintf('%s_T%03d-%03d_n%g', chave, round(1000*xtr(1)), round(1000*xtr(2)), ncr);
    end
    chave = regexprep(chave, '[^A-Za-z0-9_\-\.]', '_');
    arq_mat = fullfile(pasta, [chave '.mat']);

    if ~forcar && exist(arq_mat,'file') == 2
        S = load(arq_mat, 'Pol', 'info');
        Pol = S.Pol;  info = S.info;  info.do_estoque = true;
        if verb
            fprintf('    polar %-14s Re %.3g : ESTOQUE (%d pts, %d recuperados, %d faltando)\n', ...
                    perfil, Rer, numel(Pol.alpha), numel(info.recuperados), numel(info.faltando));
        end
        return
    end

    t0 = tic;
    if verb; fprintf('    polar %-14s Re %.3g : rodando XFOIL...\n', perfil, Rer); end
    Pol = run_xfoil(perfil, Rer, aa, struct('viscoso',visc,'n_paineis',npan, ...
                                          'ancora',anc,'dens_ba',dba,'raz_bf',rbf, ...
                                          'xtr',xtr,'ncrit',ncr));
    n_normal   = numel(Pol.alpha);
    faltando   = Pol.alpha_faltando(:).';
    recuperados = [];

    if recup && visc && ~isempty(faltando) && ~isempty(Pol.alpha)
        if verb
            fprintf('      continuacao em alpha para %d alpha(s) que nao convergiram...\n', ...
                    numel(faltando));
        end
        restantes = faltando;
        for volta = 1:numel(faltando)
            if isempty(restantes); break; end
            conv = Pol.alpha(:).';
            dist = arrayfun(@(a) min(abs(conv - a)), restantes);
            [~, idx] = min(dist);
            a = restantes(idx);
            [~, ja] = min(abs(conv - a));
            D = caminhar_alpha(arq_dat, Rer, conv(ja), a, passo_r, npan, dba, rbf, iter_r, xtr, ncr);
            novos = 0;
            for r = 1:size(D,1)
                if any(abs(Pol.alpha - D(r,1)) < 0.01); continue; end
                Pol.alpha(end+1,1) = D(r,1);  Pol.Cl(end+1,1) = D(r,2);
                Pol.Cd(end+1,1)    = D(r,3);  Pol.Cm(end+1,1) = D(r,4);
                novos = novos + 1;
                if any(abs(faltando - D(r,1)) < 0.01); recuperados(end+1) = D(r,1); end
            end
            [Pol.alpha, io] = sort(Pol.alpha);
            Pol.Cl = Pol.Cl(io);  Pol.Cd = Pol.Cd(io);  Pol.Cm = Pol.Cm(io);
            if verb
                fprintf('        %+6.2f  (de %+6.2f): %d ponto(s) novo(s)\n', a, conv(ja), novos);
            end
            restantes = setdiff(restantes, round(recuperados*100)/100);
            restantes = setdiff(restantes, a);
        end
        Pol.alpha_faltando = setdiff(round(faltando*2)/2, round(recuperados*2)/2);
    end

    suspeitos = zeros(0,3);
    normais = setdiff(round(Pol.alpha*100)/100, round(recuperados*100)/100);
    for a = recuperados
        k  = find(abs(Pol.alpha - a) < 0.01, 1);
        ab = normais(normais < a);  ac = normais(normais > a);
        if isempty(ab) || isempty(ac) || isempty(k); continue; end
        c1 = Pol.Cl(abs(Pol.alpha - ab(end)) < 0.01);  c2 = Pol.Cl(abs(Pol.alpha - ac(1)) < 0.01);
        esp = c1(1) + (c2(1) - c1(1)) * (a - ab(end)) / (ac(1) - ab(end));
        if abs(Pol.Cl(k) - esp) > 0.15
            suspeitos(end+1,:) = [a, Pol.Cl(k), esp];
        end
    end

    info = struct('perfil',perfil,'Re',Rer,'md5',hd,'chave',chave, ...
                  'n_pedidos',numel(aa),'n_normal',n_normal, ...
                  'recuperados',recuperados,'faltando',Pol.alpha_faltando(:).', ...
                  'suspeitos',suspeitos,'segundos',toc(t0), ...
                  'data',datestr(now),'do_estoque',false);
    save(arq_mat, 'Pol', 'info');

    if verb
        fprintf('      -> %d/%d pelo caminho normal, +%d por continuacao, %d ainda faltando (%.0f s)\n', ...
                n_normal, numel(aa), numel(recuperados), numel(info.faltando), info.segundos);
        for i = 1:size(suspeitos,1)
            fprintf('      [?] alpha %+5.1f recuperado com Cl %.3f; vizinhos sugerem %.3f\n', suspeitos(i,:));
        end
    end
end

function D = caminhar_alpha(arq_dat, Re, a_de, a_ate, passo, npan, dba, rbf, iter, xtr, ncr)
    D = zeros(0,4);
    tag = sprintf('%06d', randi(999999));
    cmd = ['xrec_cmd_' tag '.txt'];  pol = ['xrec_pol_' tag '.txt'];  lg = ['xrec_log_' tag '.txt'];
    apagar({cmd, pol, lg});

    n = max(1, ceil(abs(a_ate - a_de)/max(passo,1e-6)));
    seq = a_de + (a_ate - a_de)*(1:n)/n;

    fid = fopen(cmd,'w');
    fprintf(fid,'PLOP\nG F\n\n');
    fprintf(fid,'LOAD %s\n', arq_dat);
    fprintf(fid,'PPAR\nN %d\nP %.2f\nT %.2f\n\n\n', npan, dba, rbf);
    fprintf(fid,'OPER\n');
    fprintf(fid,'VISC %f\n', Re);
    if any(abs(xtr(:).' - [1 1]) > 1e-9) || abs(ncr - 9) > 1e-9
        fprintf(fid,'VPAR\nN %g\nXTR %.4f %.4f\n\n', ncr, xtr(1), xtr(2));
    end
    fprintf(fid,'ITER %d\n', iter);
    fprintf(fid,'PACC\n%s\n\n', pol);
    fprintf(fid,'ALFA %.3f\n', a_de);
    for s2 = seq; fprintf(fid,'ALFA %.3f\n', s2); end
    fprintf(fid,'PACC\n');
    fprintf(fid,'\n\n\nQUIT\n\n');
    fclose(fid);

    xfoil_exec(cmd, lg, 90);

    if exist(pol,'file') == 2
        txt = fileread(pol);
        m = regexp(txt, 'Re\s*=\s*([\d\.]+)\s*e\s*(\d+)', 'tokens', 'once');
        if ~isempty(m) && abs(str2double(m{1})*10^str2double(m{2})/Re - 1) < 0.01
            f = fopen(pol,'r');
            for i = 1:12; fgetl(f); end
            M = fscanf(f,'%f %f %f %f %f %f %f',[7,inf])';
            fclose(f);
            if ~isempty(M); D = [M(:,1), M(:,2), M(:,3), M(:,5)]; end
        end
    end
    apagar({cmd, pol, lg});
end

function h = md5_arquivo(arq)
    fid = fopen(arq,'r');  b = fread(fid, inf, '*int8');  fclose(fid);
    h = md5_bytes(b);
end

function h = md5_texto(s)
    h = md5_bytes(typecast(uint8(s(:)), 'int8'));
end

function h = md5_bytes(b)
    md = java.security.MessageDigest.getInstance('MD5');
    md.update(b);
    d = typecast(md.digest, 'uint8');
    h = lower(reshape(dec2hex(d, 2).', 1, []));
end

function apagar(lista)
    for f = lista
        if exist(f{1},'file'); delete(f{1}); end
    end
end

function v = obter(s, nome, padrao)
    if isfield(s,nome) && ~isempty(s.(nome)); v = s.(nome); else; v = padrao; end
end
