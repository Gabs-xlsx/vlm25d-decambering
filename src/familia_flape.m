function Fam = familia_flape(perfil, Re, alphas, deltas, x_dob, n_sec, opts)
%FAMILIA_FLAPE  Polares XFOIL do perfil REALMENTE defletido, uma familia por delta.
%   Fam = familia_flape('s1223.dat', Re, -8:0.5:20, [-20 -10 10 20], 0.80, 2)

    if nargin < 5 || isempty(x_dob); x_dob = 0.80;   end
    if nargin < 6 || isempty(n_sec); n_sec = 2;      end
    if nargin < 7; opts = struct(); end
    vit    = obter(opts,'viterna',true);
    oxf    = obter(opts,'xfoil',struct());
    pasta  = obter(opts,'pasta','Perfis');
    limpar = obter(opts,'limpar',false);

    if ischar(perfil) || isstring(perfil); perfil = {char(perfil)}; end
    T = teoria_flape(x_dob);
    deltas = deltas(:)';
    Fam = repmat(struct('delta',0,'Polares',[],'origem','','geom',[], ...
                        'alpha_L0',NaN,'dalphaL0_medido',NaN, ...
                        'dalphaL0_teorico',NaN,'efetividade',NaN, ...
                        'dCm_medido',NaN,'dCm_teorico',NaN, ...
                        'rms_ajuste',NaN,'n_ajuste',0,'ok',false), 1, numel(deltas));
    lixo = {};

    fprintf('  [familia] limpo ... ');  t0 = tic;
    P0 = run_xfoil(perfil{1}, Re, alphas, oxf);
    a0_limpo = alpha_L0_de(P0.alpha, P0.Cl);
    cm0_limpo = interp_seg(P0.alpha, P0.Cm, 0);
    fprintf('alpha_L0 %+.3f deg | Cm(0) %+.4f | %.0f s\n', a0_limpo, cm0_limpo, toc(t0));

    for m = 1:numel(deltas)
        d = deltas(m);
        Fam(m).delta = d;
        Fam(m).dalphaL0_teorico = T.dalphaL0(d);
        Fam(m).dCm_teorico      = T.m2 * tan(deg2rad(d));

        if abs(d) < 1e-9
            Pol = P0;
        else
            tag = sprintf('_f%+03.0f', d);  tag = strrep(tag,'+','p');  tag = strrep(tag,'-','m');
            [~,b,e] = fileparts(perfil{1});
            arq = [b tag e];
            G = flapear_perfil(perfil{1}, d, x_dob, arq, ...
                    struct('pasta',pasta,'raio',obter(opts,'raio',0),'verbose',true));
            Fam(m).geom = rmfield(G, {'XY','XY_limpo','XY_original'});
            lixo{end+1} = fullfile(pasta, arq);
            if ~G.ok
                warning('familia_flape:geom','geometria de delta %+.1f suspeita -- seguindo mesmo assim.', d);
            end
            sh  = Fam(m).dalphaL0_teorico;
            passo = median(diff(sort(alphas(:))));
            if ~isfinite(passo) || passo <= 0; passo = 0.5; end
            piso = obter(opts,'alpha_piso',-16);
            teto = obter(opts,'alpha_teto', 18);
            lo  = max(piso, min(alphas) + sh);
            hi  = min(teto, max(alphas) + min(sh,0));
            aa  = lo:passo:hi;
            oxd = oxf;
            anc_base = 4;
            if isfield(oxf,'ancora') && ~isempty(oxf.ancora); anc_base = oxf.ancora; end
            if isfield(oxf,'ancora') && isempty(oxf.ancora)
                txt_anc = 'sem (crescente pura)';
            else
                oxd.ancora = min(max(anc_base + sh, lo + 1), min(hi - 1, 8));
                txt_anc = sprintf('%+.1f', oxd.ancora);
            end
            fprintf('  [familia] delta %+6.1f | XFOIL de %+.1f a %+.1f, ancora %s ... ', ...
                    d, min(aa), max(aa), txt_anc);  t0 = tic;
            Pol = run_xfoil(arq, Re, aa, oxd);
            fprintf('%d pts | %.0f s\n', numel(Pol.alpha), toc(t0));
        end

        if isempty(Pol.alpha)
            warning('familia_flape:vazia','delta %+.1f: polar vazia.', d);
            Fam(m).Polares = [];  continue;
        end

        [dsh, dcm, rmsq, npt] = deslocamento_polar(P0.alpha, P0.Cl, P0.Cm, ...
                                                   Pol.alpha, Pol.Cl, Pol.Cm);
        Fam(m).dalphaL0_medido = dsh;
        Fam(m).dCm_medido      = dcm;
        Fam(m).rms_ajuste      = rmsq;
        Fam(m).n_ajuste        = npt;
        Fam(m).alpha_L0        = a0_limpo + dsh;
        if abs(d) > 1e-9
            Fam(m).efetividade = Fam(m).dalphaL0_medido / Fam(m).dalphaL0_teorico;
        else
            Fam(m).efetividade = 1;
        end

        if vit
            Pv = extrapolar_viterna(Pol);
        else
            Pv = Pol;
        end
        Pv.flape = struct('delta_deg',d,'x_dob',x_dob,'fonte','XFOIL');
        Fam(m).Polares = repmat(Pv, 1, n_sec);
        Fam(m).origem  = 'XFOIL defletido';

        bom_ef = Fam(m).efetividade > 0.15 && Fam(m).efetividade < 1.30;
        bom_n  = npt >= 20 && numel(Pol.alpha) >= 30;
        Fam(m).ok = abs(d) < 1e-9 || (bom_ef && bom_n);
        fprintf(['             alpha_L0 %+.3f (limpo %+.3f) | d_medido %+.3f, ' ...
                 'd_teorico %+.3f -> efetividade %.2f | dCm %+.4f vs %+.4f | ' ...
                 'ajuste rms %.4f em %d pts | %s\n'], ...
                Fam(m).alpha_L0, a0_limpo, Fam(m).dalphaL0_medido, ...
                Fam(m).dalphaL0_teorico, Fam(m).efetividade, ...
                Fam(m).dCm_medido, Fam(m).dCm_teorico, rmsq, npt, ...
                char(79*Fam(m).ok + 88*~Fam(m).ok));
    end

    for s = [-1 1]
        idx = find(sign([Fam.delta]) == s);
        if numel(idx) < 2; continue; end
        [~, o] = sort(abs([Fam(idx).delta]));
        idx = idx(o);
        for q = 2:numel(idx)
            a = idx(q-1);  b = idx(q);
            dm = abs(Fam(b).dalphaL0_medido) - abs(Fam(a).dalphaL0_medido);
            dt = abs(Fam(b).dalphaL0_teorico) - abs(Fam(a).dalphaL0_teorico);
            if ~isfinite(dm) || dt <= 0; continue; end
            if dm < 0.25*dt
                Fam(b).ok = false;
                warning('familia_flape:monotonia', ...
                  ['delta %+.1f: o deslocamento medido cresceu so %.2f deg sobre ' ...
                   'delta %+.1f, contra %.2f deg de teoria (%.0f%%). Familia ' ...
                   'reprovada -- dado insuficiente, nao efetividade baixa.'], ...
                  Fam(b).delta, dm, Fam(a).delta, dt, 100*dm/dt);
            end
        end
    end

    if limpar
        for i = 1:numel(lixo)
            if exist(lixo{i},'file'); delete(lixo{i}); end
        end
    end
end

function [dsh, dcm, rmsq, npt] = deslocamento_polar(ac, clc, cmc, af, clf, cmf)

    dsh = NaN;  dcm = NaN;  rmsq = NaN;  npt = 0;
    ac = ac(:); clc = clc(:); cmc = cmc(:);
    af = af(:); clf = clf(:); cmf = cmf(:);
    if numel(ac) < 8 || numel(af) < 8; return; end

    [~,ic] = max(clc);  [~,ifl] = max(clf);
    mc = ac <= ac(ic) - 1;   mf = af <= af(ifl) - 1;
    if sum(mc) < 8 || sum(mf) < 8; return; end
    ac = ac(mc); clc = clc(mc); cmc = cmc(mc);
    af = af(mf); clf = clf(mf); cmf = cmf(mf);

    minimo = max(8, ceil(0.5*numel(af)));
    D = -25:0.05:25;
    e = nan(size(D));  n = zeros(size(D));
    for i = 1:numel(D)
        aq = af - D(i);
        k  = aq >= min(ac) & aq <= max(ac);
        if sum(k) < minimo; continue; end
        r = clf(k) - interp1(ac, clc, aq(k), 'linear');
        e(i) = sqrt(mean(r.^2));  n(i) = sum(k);
    end
    if all(isnan(e)); return; end

    [~,i0] = min(e);
    dsh = D(i0);  rmsq = e(i0);  npt = n(i0);
    if i0 > 1 && i0 < numel(D) && all(isfinite(e([i0-1 i0 i0+1])))
        y1 = e(i0-1); y2 = e(i0); y3 = e(i0+1);
        den = y1 - 2*y2 + y3;
        if abs(den) > eps
            dsh = D(i0) + 0.5*(y1 - y3)/den * (D(2)-D(1));
        end
    end

    aq = af - dsh;
    k  = aq >= min(ac) & aq <= max(ac);
    if sum(k) >= minimo
        dcm = mean(cmf(k) - interp1(ac, cmc, aq(k), 'linear'));
    end
end

function a0 = alpha_L0_de(a, cl)
    a0 = NaN;  a = a(:);  cl = cl(:);
    s = sign(cl);
    k = find(s(1:end-1) < 0 & s(2:end) >= 0);
    if isempty(k); return; end
    [~,q] = min(abs(a(k)));  i = k(q);
    if cl(i+1) == cl(i); a0 = a(i); return; end
    a0 = a(i) - cl(i)*(a(i+1)-a(i))/(cl(i+1)-cl(i));
end

function v = interp_seg(x, y, xq)
    if numel(x) < 2; v = NaN; return; end
    v = interp1(x, y, xq, 'linear', 'extrap');
end

function v = obter(s,nome,padrao)
    if isstruct(s) && isfield(s,nome) && ~isempty(s.(nome)); v = s.(nome); else; v = padrao; end
end
