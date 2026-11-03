function [xc, zc] = perfil_linha_media(nome_perfil, n_pontos)
% PERFIL_LINHA_MEDIA  Linha media (camber) adimensional de um perfil .dat.

    if nargin < 2 || isempty(n_pontos); n_pontos = 101; end

    persistent cache
    if isempty(cache); cache = struct('chave', {}, 'xc', {}, 'zc', {}); end
    chave = sprintf('%s|%d', lower(nome_perfil), n_pontos);
    for k = 1:numel(cache)
        if strcmp(cache(k).chave, chave)
            xc = cache(k).xc; zc = cache(k).zc; return;
        end
    end

    caminho = resolver_caminho_perfil(nome_perfil);

    fid = fopen(caminho, 'r');
    if fid < 0; error('Falha ao abrir %s', caminho); end
    limpar = onCleanup(@() fclose(fid));
    primeira = fgetl(fid);
    if ~isempty(str2num(primeira))
        frewind(fid);
    end
    dados = fscanf(fid, '%f %f', [2, inf])';
    if size(dados,1) < 10
        error('Arquivo de perfil %s tem poucos pontos (%d).', caminho, size(dados,1));
    end

    x = dados(:,1);  z = dados(:,2);
    [~, i_ba] = min(x);

    x_sup = x(1:i_ba);       z_sup = z(1:i_ba);
    x_inf = x(i_ba:end);     z_inf = z(i_ba:end);
    [x_sup, i] = unique(x_sup); z_sup = z_sup(i);
    [x_inf, i] = unique(x_inf); z_inf = z_inf(i);

    xc = 0.5 * (1 - cos(linspace(0, pi, n_pontos)'));
    zc = ( interp1(x_sup, z_sup, xc, 'pchip', 'extrap') + ...
           interp1(x_inf, z_inf, xc, 'pchip', 'extrap') ) / 2;

    cache(end+1) = struct('chave', chave, 'xc', xc, 'zc', zc);
end
