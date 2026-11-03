function R = flapear_perfil(arquivo, delta_deg, x_dob, saida, opts)
%FLAPEAR_PERFIL  Gira o bordo de fuga de um perfil 2D e grava o .dat defletido.
%   R = flapear_perfil('s1223.dat', 10, 0.80, 's1223_f10.dat')

    if nargin < 5; opts = struct(); end
    pasta = obter(opts,'pasta','Perfis');
    raio  = obter(opts,'raio',0);
    verb  = obter(opts,'verbose',true);

    ent = fullfile(pasta, arquivo);
    if exist(ent,'file') ~= 2
        error('flapear_perfil:perfil','Perfil "%s" nao encontrado.', ent);
    end

    [XY, nome] = ler_dat(ent);
    if size(XY,1) < 20
        error('flapear_perfil:poucos','So li %d pontos de "%s" -- formato nao reconhecido.', ...
              size(XY,1), ent);
    end
    XY0 = XY;

    XY = para_selig(XY);

    [~, iLE] = min(XY(:,1));
    sup = XY(1:iLE,   :);
    inf_= XY(iLE:end, :);

    yh_s = interp1(sup(:,1), sup(:,2), x_dob, 'linear');
    yh_i = interp1(inf_(:,1), inf_(:,2), x_dob, 'linear');
    if ~isfinite(yh_s) || ~isfinite(yh_i)
        error('flapear_perfil:fora','x_dob = %.3f cai fora das coordenadas do perfil.', x_dob);
    end
    if min(abs(sup(:,1)-x_dob)) > 1e-7
        k = find(sup(:,1) < x_dob, 1, 'first');
        sup = [sup(1:k-1,:); x_dob yh_s; sup(k:end,:)];
    end
    if min(abs(inf_(:,1)-x_dob)) > 1e-7
        k = find(inf_(:,1) > x_dob, 1, 'first');
        inf_ = [inf_(1:k-1,:); x_dob yh_i; inf_(k:end,:)];
    end
    XY = [sup; inf_(2:end,:)];

    y_dob = 0.5*(yh_s + yh_i);

    n_gir = 0;
    XYd = XY;
    girado = false(size(XY,1),1);
    for i = 1:size(XY,1)
        dx = XY(i,1) - x_dob;
        if dx <= 1e-12; continue; end
        girado(i) = true;
        if raio > 0
            s = min(dx/raio, 1);  w = s*s*(3 - 2*s);
        else
            w = 1;
        end
        ang = -deg2rad(delta_deg) * w;
        c = cos(ang);  s2 = sin(ang);
        v = [dx, XY(i,2) - y_dob];
        XYd(i,:) = [x_dob, y_dob] + [c*v(1) - s2*v(2), s2*v(1) + c*v(2)];
        n_gir = n_gir + 1;
    end

    invade = girado & (XYd(:,1) < x_dob - 1e-9);
    n_sel = sum(invade);
    if n_sel > 0
        XYd(invade,:) = [];
        XY(invade,:)  = [];
        girado(invade) = [];
        n_gir = n_gir - n_sel;
    end

    R = struct();
    R.arquivo_entrada = ent;
    R.delta_deg = delta_deg;
    R.x_dob = x_dob;
    R.y_dob = y_dob;
    R.n_pontos = size(XYd,1);
    R.n_selados = n_sel;
    R.n_girados = n_gir;

    ibf0 = 1;
    R.dz_bf = XYd(ibf0,2) - XY(ibf0,2);
    R.dz_bf_esperado = -sin(deg2rad(delta_deg)) * (XY(ibf0,1) - x_dob) ...
                       + (cos(deg2rad(delta_deg)) - 1) * (XY(ibf0,2) - y_dob);
    R.erro_bf = abs(R.dz_bf - R.dz_bf_esperado);
    R.dz_bf_simples = -(1 - x_dob) * sin(deg2rad(delta_deg));

    R.corda_antes  = max(XY(:,1))  - min(XY(:,1));
    R.corda_depois = max(XYd(:,1)) - min(XYd(:,1));
    R.esp_max_antes  = espessura_max(XY);
    R.esp_max_depois = espessura_max(XYd);

    T = teoria_flape(x_dob);
    R.c2 = T.c2;  R.m2 = T.m2;
    R.dalphaL0_teorico = T.dalphaL0(delta_deg);
    R.dCm_teorico      = T.m2 * tan(deg2rad(delta_deg));

    if nargin >= 4 && ~isempty(saida)
        sai = fullfile(pasta, saida);
        fid = fopen(sai,'w');
        if fid < 0; error('flapear_perfil:grava','nao consegui escrever "%s".', sai); end
        fprintf(fid, '%s flape %+.1f deg em x/c %.3f\n', nome, delta_deg, x_dob);
        fprintf(fid, '%12.7f %12.7f\n', XYd');
        fclose(fid);
        R.arquivo_saida = sai;
        XYv = ler_dat(sai);
        R.erro_ida_volta = max(max(abs(XYv - XYd)));
    else
        R.arquivo_saida = '';
        R.erro_ida_volta = 0;
    end
    R.XY = XYd;
    R.XY_limpo = XY;
    R.XY_original = XY0;

    R.ok = R.erro_bf < 1e-9 && R.erro_ida_volta < 1e-6 && ...
           abs(R.esp_max_depois - R.esp_max_antes) < 5e-3 && ...
           (abs(delta_deg) < 1e-9 || n_gir > 5);

    if verb
        fprintf('flape %+6.1f deg | girados %3d/%3d pts | dz_BF %+.5f (esperado %+.5f, forma fechada %+.5f, erro %.1e)\n', ...
                delta_deg, n_gir, R.n_pontos, R.dz_bf, R.dz_bf_esperado, ...
                R.dz_bf_simples, R.erro_bf);
        fprintf('             | esp_max %.5f -> %.5f | ida-volta %.1e | teoria: d_alphaL0 %+.3f deg, d_Cm %+.4f | %s\n', ...
                R.esp_max_antes, R.esp_max_depois, R.erro_ida_volta, ...
                R.dalphaL0_teorico, R.dCm_teorico, char(79*R.ok + 88*~R.ok));
    end
    if ~R.ok
        warning('flapear_perfil:suspeito', ...
            'geometria defletida NAO passou na conferencia (delta = %+.1f).', delta_deg);
    end
end

function [XY, nome] = ler_dat(caminho)
    txt = fileread(caminho);
    L = regexp(txt, '\r\n|\n|\r', 'split');
    XY = zeros(0,2);  nome = '';
    for i = 1:numel(L)
        s = strtrim(L{i});
        if isempty(s); continue; end
        v = sscanf(s, '%f %f');
        if numel(v) == 2 && all(isfinite(v)) && abs(v(1)) <= 2 && abs(v(2)) <= 2
            XY(end+1,:) = v(:)';
        elseif isempty(nome)
            nome = s;
        end
    end
    if isempty(nome); nome = 'perfil'; end
end

function XY = para_selig(XY)
    if XY(1,1) > 0.5; return; end
    d = diff(XY(:,1));
    q = find(d < -1e-9, 1, 'first');
    if isempty(q)
        XY = flipud(XY);  return;
    end
    sup = XY(1:q,   :);
    inf_= XY(q+1:end,:);
    XY  = [flipud(sup); inf_];
end

function e = espessura_max(XY)
    [~,iLE] = min(XY(:,1));
    s = flipud(XY(1:iLE,:));
    t = XY(iLE:end,:);
    [xs, is] = unique(s(:,1), 'stable');  ys = s(is,2);
    [xt, it] = unique(t(:,1), 'stable');  yt = t(it,2);
    if numel(xs) < 3 || numel(xt) < 3; e = NaN; return; end
    lo = max(min(xs), min(xt));  hi = min(max(xs), max(xt));
    xg = linspace(lo + 1e-6, hi - 1e-6, 200);
    e = max(interp1(xs, ys, xg, 'linear') - interp1(xt, yt, xg, 'linear'));
end

function v = obter(s,nome,padrao)
    if isfield(s,nome) && ~isempty(s.(nome)); v = s.(nome); else; v = padrao; end
end
