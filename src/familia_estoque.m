function FAMX = familia_estoque(perfil, Re, deltas, x_dob, n_sec, xtr, ncrit, forcar, limite_s, alfas)
%FAMILIA_ESTOQUE  Familia de polares do perfil REALMENTE defletido, com cache.
%   FAMX = familia_estoque('s1223.dat', 250e3, [-15 -10 -5 5 10 15], 0.75, 2)
%   FAMX = familia_estoque(..., [0.02 0.02], 9)      % transicao forcada

    if nargin < 6 || isempty(xtr);   xtr = [1 1]; end
    if nargin < 7 || isempty(ncrit); ncrit = 9;   end
    if nargin < 8 || isempty(forcar); forcar = false; end
    if nargin < 9 || isempty(limite_s); limite_s = 900; end
    if nargin < 10 || isempty(alfas); alfas = -12:0.5:16; end

    [~, raiz] = fileparts(perfil);
    tag = '';
    if any(abs(xtr(:).' - [1 1]) > 1e-9)
        tag = sprintf('_T%03d-%03d', round(1000*xtr(1)), round(1000*xtr(2)));
    end
    if abs(ncrit - 9) > 1e-9; tag = [tag sprintf('_n%g', ncrit)]; end
    if exist('Resultados_VLM','dir') ~= 7; mkdir('Resultados_VLM'); end
    arq = fullfile('Resultados_VLM', sprintf('familia_%s_Re%.0f%s_x%03d.mat', ...
                   raiz, Re, tag, round(1000*x_dob)));

    deltas = unique(deltas(:)');
    deltas = deltas(abs(deltas) > 1e-9);
    if isempty(deltas); FAMX = []; return; end

    FAMX = [];
    if ~forcar && exist(arq,'file') == 2
        Z = load(arq,'FAMX');
        if isfield(Z,'FAMX'); FAMX = Z.FAMX; end
    end

    if isempty(FAMX)
        ja = [];
    else
        cheias = arrayfun(@(f) ~isempty(f.Polares), FAMX);
        ja = [FAMX(cheias).delta];
        if any(~cheias)
            fprintf('    familia: %s em cache como vazia, sera refeita\n', ...
                    mat2str([FAMX(~cheias).delta]));
            FAMX = FAMX(cheias);
        end
    end
    falta = deltas(~ismember(round(deltas*1e3), round(ja*1e3)));

    if ~isempty(falta)
        fprintf('    familia XFOIL %s Re %.0f%s: rodando %s na geometria defletida...\n', ...
                raiz, Re, tag, mat2str(falta));
        t0 = tic;
        Fn = familia_flape(perfil, Re, alfas, falta, x_dob, n_sec, ...
                           struct('xfoil', struct('xtr',xtr,'ncrit',ncrit,'limite_s',limite_s)));
        if isempty(FAMX); FAMX = Fn(:)'; else; FAMX = unir_fam(FAMX(:)', Fn(:)'); end
        save(arq, 'FAMX', '-v7.3');
        fprintf('    familia XFOIL: %.1f min -> %s\n', toc(t0)/60, arq);
    end

    [~, iq] = sort([FAMX.delta]);  FAMX = FAMX(iq);
    FAMX = FAMX(ismember(round([FAMX.delta]*1e3), round(deltas*1e3)));
    for k = 1:numel(FAMX)
        if isempty(FAMX(k).Polares)
            warning('familia_estoque:vazia', ...
                    'familia vazia em delta %+g (Re %.0f%s): o XFOIL nao fechou a geometria defletida.', ...
                    FAMX(k).delta, Re, tag);
        end
    end
end

function C = unir_fam(A, B)
    ca = fieldnames(A);  cb = fieldnames(B);
    for k = 1:numel(cb)
        if ~isfield(A, cb{k}); [A.(cb{k})] = deal([]); end
    end
    for k = 1:numel(ca)
        if ~isfield(B, ca{k}); [B.(ca{k})] = deal([]); end
    end
    B = orderfields(B, A);
    C = [A(:)', B(:)'];
end
