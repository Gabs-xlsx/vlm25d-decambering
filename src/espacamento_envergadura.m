function Y_nodes = espacamento_envergadura(quebras, Ny_alvo, tipo, aperto)
% ESPACAMENTO_ENVERGADURA  Nos de envergadura alinhados com as quebras, com a

    if nargin < 3 || isempty(tipo);   tipo = 1;      end
    if nargin < 4 || isempty(aperto); aperto = 0.75; end
    aperto = min(max(aperto, 0), 1);
    if tipo == 0; aperto = 0; end

    quebras = quebras(:)';
    y0 = quebras(1);
    L  = quebras(end) - quebras(1);
    if L <= 0
        error('espacamento_envergadura:quebras', 'quebras degeneradas.');
    end

    if tipo == 2
        Y_nodes = quebras(1);
        for k = 1:(numel(quebras)-1)
            span_k = quebras(k+1) - quebras(k);
            n_k    = max(1, round(Ny_alvo * span_k / L));
            t      = distribuir_espacamento(n_k, 1);
            Y_nodes = [Y_nodes, quebras(k) + span_k * t(2:end)];
        end
        return
    end

    s_q = mapa_inverso((quebras - y0)/L, aperto);
    s_q(1) = 0;  s_q(end) = 1;

    Y_nodes = quebras(1);
    for k = 1:(numel(quebras)-1)
        ds  = s_q(k+1) - s_q(k);
        n_k = max(1, round(Ny_alvo * ds));
        s_k = s_q(k) + ds * (1:n_k) / n_k;
        y_k = y0 + L * mapa_direto(s_k, aperto);
        y_k(end) = quebras(k+1);
        Y_nodes = [Y_nodes, y_k];
    end
end

function u = mapa_direto(s, w)
    u = (1-w)*s + w*(1 - cos(pi*s))/2;
end

function s = mapa_inverso(u, w)
    u = min(max(u(:)', 0), 1);
    s = u;
    for it = 1:60
        f  = mapa_direto(s, w) - u;
        d  = (1-w) + w*(pi/2)*sin(pi*s);
        s  = min(max(s - f./max(d, 1e-9), 0), 1);
    end
end
