function caminho = resolver_caminho_perfil(nome_perfil)
% RESOLVER_CAMINHO_PERFIL  Fonte unica de verdade para localizar um .dat de perfil.

    if exist(nome_perfil, 'file') == 2
        caminho = nome_perfil; return;
    end

    [~, base, ext] = fileparts(nome_perfil);
    alvo = lower([base ext]);

    for pasta = {'Perfis', '.'}
        p = pasta{1};
        if exist(p, 'dir') ~= 7; continue; end
        lista = dir(fullfile(p, '*.dat'));
        for k = 1:numel(lista)
            if strcmp(lower(lista(k).name), alvo)
                caminho = fullfile(p, lista(k).name); return;
            end
        end
    end

    error('resolver_caminho_perfil:naoEncontrado', ...
          ['Perfil "%s" nao encontrado em Perfis/ nem na raiz do projeto.\n' ...
           'Diretorio atual: %s'], nome_perfil, pwd);
end
