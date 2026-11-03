function [S, D] = mapa_liberacao(resolvidas, dj, yf, modo, t_max, n_aj)
%MAPA_LIBERACAO  Matriz que espalha o estado de decambering das faixas

    resolvidas = logical(resolvidas(:));
    dj = dj(:);  yf = yf(:);
    Ny = numel(resolvidas);
    if nargin < 4 || isempty(modo);  modo  = 'copia'; end
    if nargin < 5 || isempty(t_max); t_max = 2;       end
    if nargin < 6 || isempty(n_aj);  n_aj  = 4;       end

    idx_res = find(resolvidas);  M = numel(idx_res);
    if M == 0
        error('mapa_liberacao:sem_incognita', ...
              'as %d faixas foram liberadas -- nao sobrou incognita.', Ny);
    end
    pos = zeros(Ny,1);  pos(idx_res) = 1:M;

    seg = ones(Ny,1);
    for j = 2:Ny
        seg(j) = seg(j-1) + double(abs(dj(j) - dj(j-1)) > 1e-6);
    end
    ns = seg(end);
    n_res_seg = accumarray(seg, double(resolvidas), [ns 1]);
    seg_vazio = find(n_res_seg == 0).';

    seg_n   = accumarray(seg, 1, [ns 1]);
    seg_lib = accumarray(seg, double(~resolvidas), [ns 1]);
    D = struct('modo',modo,'seg',seg,'n_seg',ns,'n_res_seg',n_res_seg(:).', ...
               'seg_vazio',seg_vazio,'origem',zeros(Ny,1),'t',zeros(Ny,1), ...
               'cruzou',0,'t_max',t_max,'n_ajuste',n_aj, ...
               'seg_n',seg_n(:).', 'seg_lib',seg_lib(:).', ...
               'seg_frac',(seg_lib./max(seg_n,1)).', 'frac_max',max(seg_lib./max(seg_n,1)));

    lin = zeros((n_aj+1)*Ny,1);  col = lin;  val = lin;  nz = 0;

    switch lower(modo)
    case 'copia'
        for j = 1:Ny
            [~,q] = min(abs(idx_res - j));
            nz = nz + 1;  lin(nz) = j;  col(nz) = q;  val(nz) = 1;
            if ~resolvidas(j)
                D.origem(j) = 1;
                if seg(idx_res(q)) ~= seg(j); D.cruzou = D.cruzou + 1; end
            end
        end
    case 'segmento'
        if ~isempty(seg_vazio)
            error('mapa_liberacao:segmento_vazio', ...
                  ['segmento(s) %s sem nenhuma faixa resolvida: uma superficie ' ...
                   'inteira foi liberada e herdaria o estado de outra deflexao. ' ...
                   'Reduza n_borda/n_extrap ou refine a malha.'], mat2str(seg_vazio));
        end
        for j = 1:Ny
            if resolvidas(j)
                nz = nz + 1;  lin(nz) = j;  col(nz) = pos(j);  val(nz) = 1;
                continue
            end
            rs = idx_res(seg(idx_res) == seg(j));
            ab = rs(rs < j);   ac = rs(rs > j);

            if ~isempty(ab) && ~isempty(ac)
                a = ab(end);  b = ac(1);  w = 0.5;
                if abs(yf(b) - yf(a)) > 1e-12
                    w = (yf(j) - yf(a)) / (yf(b) - yf(a));
                end
                w = min(max(w,0),1);
                nz=nz+1; lin(nz)=j; col(nz)=pos(a); val(nz)=1-w;
                nz=nz+1; lin(nz)=j; col(nz)=pos(b); val(nz)=w;
                D.origem(j) = 2;

            elseif numel(rs) >= 2
                if isempty(ac); ord = rs(end:-1:1); else; ord = rs(:); end
                P = ord(1:min(n_aj, numel(ord)));
                k = numel(P);
                yp = yf(P);  ybar = mean(yp);  Sxx = sum((yp - ybar).^2);
                if Sxx < 1e-20
                    nz=nz+1; lin(nz)=j; col(nz)=pos(P(1)); val(nz)=1;
                    D.origem(j) = 1;
                else
                    bs = yf(P(1)) - yf(P(end));  t = 0;
                    if abs(bs) > 1e-12; t = (yf(j) - yf(P(1))) / bs; end
                    t = min(max(t,0), t_max);
                    y_ef = yf(P(1)) + t*bs;
                    for m = 1:k
                        w = 1/k + (y_ef - ybar)*(yp(m) - ybar)/Sxx;
                        nz=nz+1; lin(nz)=j; col(nz)=pos(P(m)); val(nz)=w;
                    end
                    D.origem(j) = 3;  D.t(j) = t;
                end

            else
                nz=nz+1; lin(nz)=j; col(nz)=pos(rs(1)); val(nz)=1;
                D.origem(j) = 1;
            end
        end
    otherwise
        error('mapa_liberacao:modo', ...
              'modo "%s" desconhecido (use copia ou segmento).', modo);
    end

    S = sparse(lin(1:nz), col(1:nz), val(1:nz), Ny, M);
end
