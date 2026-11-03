function Aeronave = montar_aeronave(G)
% MONTAR_AERONAVE  Geometria parametrica no MESMO formato do gera_geometria.m.

    b2  = G.semi_env;
    cr  = G.corda_raiz;
    ct_ = obter(G,'corda_ponta',cr);
    swp = obter(G,'enflech_ba',0);
    die = obter(G,'diedro',0);
    tr  = obter(G,'torcao_raiz',0);
    tt  = obter(G,'torcao_ponta',0);
    per = obter(G,'perfil','s1223.dat');

    eta = 0;  cor = cr;
    Q = obter(G,'quebras',[]);
    for k = 1:numel(Q)
        eta(end+1) = Q(k).eta;    cor(end+1) = Q(k).corda;
    end
    eta(end+1) = 1;  cor(end+1) = ct_;
    [eta, io] = sort(eta);  cor = cor(io);

    for k = 1:numel(eta)
        y = eta(k)*b2;
        S(k) = struct('Y', y, ...
                      'X_LE', y*tand(swp), ...
                      'Z_LE', y*tand(die), ...
                      'Corda', cor(k), ...
                      'Torcao', tr + (tt - tr)*eta(k), ...
                      'Perfil', per);
    end

    Aeronave = struct();
    Aeronave.nome = obter(G,'nome','Asa_parametrica');
    Aeronave.Xcg  = obter(G,'Xcg',0.25*cr);
    Aeronave.Zcg  = obter(G,'Zcg',0);
    Aeronave.Malha = struct('nx_corda',32,'ny_env',30,'espac_x',1,'espac_y',1, ...
                            'aperto_ponta',0.75);
    Aeronave.Asa.Secoes = S;

    SF = obter(G,'superficies',[]);
    if isempty(SF)
        A = obter(G,'aileron',[]);
        if ~isempty(A)
            SF = struct('nome',obter(A,'nome','aileron'),'tipo','assimetrico', ...
                        'inicio_eta',A.inicio_eta,'fim_eta',A.fim_eta,'corda_c',A.corda_c);
        end
    end
    for k = 1:numel(SF)
        a = SF(k);
        if a.inicio_eta >= a.fim_eta || a.fim_eta > 1 || a.inicio_eta < 0
            error('montar_aeronave:superficie', ...
                  '"%s" fora da asa: inicio_eta %.3f, fim_eta %.3f', a.nome, a.inicio_eta, a.fim_eta);
        end
        for m = 1:k-1
            b = SF(m);
            if a.inicio_eta < b.fim_eta - 1e-9 && b.inicio_eta < a.fim_eta - 1e-9
                error('montar_aeronave:sobreposicao', ...
                      '"%s" e "%s" se sobrepoem na envergadura.', a.nome, b.nome);
            end
        end
        tipo = 'assimetrico';
        if isfield(a,'tipo') && ~isempty(a.tipo); tipo = a.tipo; end
        Aeronave.Asa.Controles(k) = struct('nome', a.nome, 'tipo', tipo, ...
            'inicio_y', a.inicio_eta*b2, 'fim_y', a.fim_eta*b2, 'corda_c', a.corda_c);
    end

    Aeronave.CasosVoo(1) = struct('nome','main','V_inf',12,'rho',1.225,'alpha',0);

    ref = referencias_da_geometria(Aeronave);
    Aeronave.Sref = ref.Sref;  Aeronave.Bref = ref.Bref;  Aeronave.Cref = ref.Cref;
    Aeronave.Ref  = ref;
end

function v = obter(s, nome, padrao)
    if isstruct(s) && isfield(s,nome) && ~isempty(s.(nome)); v = s.(nome); else; v = padrao; end
end
