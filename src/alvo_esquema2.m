function [Cl_t, Cm_t, a_t, estol, n_raiz] = alvo_esquema2(Pol_a, Pol_Cl, Pol_Cm, ...
                                          a_s, Cl_s, a_p, Cl_p, estol_ant, modo)
% ALVO_ESQUEMA2  Alvo viscoso pela RETA DE TRAJETORIA (M&G 2006, Esquema 2).

    if nargin < 9 || isempty(modo); modo = 'max'; end
    [Cl_max, i_max] = max(Pol_Cl);
    a_stall = Pol_a(i_max);
    n_raiz = 0;

    if abs(a_p - a_s) < 1e-9
        a_t = a_s;
        Cl_t = interp1(Pol_a, Pol_Cl, a_t, 'linear', 'extrap');
        Cm_t = interp1(Pol_a, Pol_Cm, a_t, 'linear', 'extrap');
        estol = a_t > a_stall;  n_raiz = 0;
        return;
    end

    m = (Cl_p - Cl_s) / (a_p - a_s);
    f = @(a) (Cl_s + m*(a - a_s)) - interp1(Pol_a, Pol_Cl, a, 'linear', 'extrap');

    a_scan = Pol_a(:);
    fv = arrayfun(f, a_scan);
    troca = find(fv(1:end-1) .* fv(2:end) <= 0);

    raizes = [];
    for k = troca(:)'
        aa = a_scan(k); bb = a_scan(k+1);
        if fv(k) == 0; raizes(end+1) = aa; continue; end
        for it = 1:40
            cc = 0.5*(aa+bb);
            if f(aa)*f(cc) <= 0; bb = cc; else; aa = cc; end
        end
        raizes(end+1) = 0.5*(aa+bb);
    end

    if estol_ant && a_s < a_stall - 2
        estol_ant = false;
    end

    n_raiz = numel(raizes);
    if isempty(raizes)
        a_t = a_s;
    elseif numel(raizes) == 1
        a_t = raizes(1);
    else
        pre  = raizes(raizes <= a_stall);
        pos  = raizes(raizes >  a_stall);
        if estol_ant && ~isempty(pos)
            a_t = escolher(pos, a_s, modo);
        elseif ~estol_ant && ~isempty(pre)
            a_t = escolher(pre, a_s, modo);
        else
            [~,i] = min(abs(raizes - a_s));
            a_t = raizes(i);
        end
    end

    Cl_t = interp1(Pol_a, Pol_Cl, a_t, 'linear', 'extrap');
    Cm_t = interp1(Pol_a, Pol_Cm, a_t, 'linear', 'extrap');
    estol = a_t > a_stall;
end

function a = escolher(cand, a_s, modo)
    if strcmpi(modo,'proxima')
        [~,i] = min(abs(cand - a_s));  a = cand(i);
    else
        a = max(cand);
    end
end
