function H = momento_dobradica(Malha, Res, Aeronave, nome, V, rho)
% MOMENTO_DOBRADICA  Momento de dobradica de uma superficie de comando.

    if nargin < 5 || isempty(V);   V   = 12;    end
    if nargin < 6 || isempty(rho); rho = 1.225; end
    q = 0.5*rho*V^2;

    Nx = Malha.Nx;  Ny = Malha.Ny;  nf = Malha.Nx_front;
    P  = Res.Paineis;
    if ~isfield(P,'Forca')
        error('momento_dobradica:sem_forca', ...
              'Res.Paineis nao tem o campo Forca -- passe o Res do resolver_vlm.');
    end

    H = struct('nome',nome, ...
        'Ch_dir',NaN,'Ch_esq',NaN,'Ch_dif',NaN,'Ch_tot',NaN, ...
        'H_dir',0,'H_esq',0,'Sf_dir',0,'Sf_esq',0,'cf_dir',NaN,'cf_esq',NaN, ...
        'Ch_faixa',nan(Ny,1),'H_faixa',nan(Ny,1),'Sf_faixa',nan(Ny,1),'y',nan(Ny,1), ...
        'cf_faixa',nan(Ny,1));
    Hd = 0; Sd = 0; ed = 0;
    He = 0; Se = 0; ee = 0;

    for j = 1:Ny
        base = (j-1)*Nx;
        alvo = false;
        for i = nf+1:Nx
            if P(base+i).is_control && strcmpi(P(base+i).nome_controle, nome)
                alvo = true; break;
            end
        end
        if ~alvo; continue; end

        pf = P(base + nf + 1);
        A  = pf.P1;  k = pf.P4 - pf.P1;  nk = norm(k);
        if nk < 1e-12; continue; end
        k = k / nk;

        Hj = 0;  Sj = 0;
        for i = nf+1:Nx
            p = P(base+i);
            A_no = p.P1 + 0.25*(p.P2 - p.P1);
            B_no = p.P4 + 0.25*(p.P3 - p.P4);
            r    = (A_no + B_no)/2 - A;
            Hj   = Hj + dot(cross(r, p.Forca), k);
            Sj   = Sj + p.Area;
        end
        yj = pf.CP(2);
        H.y(j) = yj/(Aeronave.Bref/2);
        H.H_faixa(j) = Hj;  H.Sf_faixa(j) = Sj;  H.cf_faixa(j) = Sj/nk;
        if yj >= 0;  Hd = Hd + Hj;  Sd = Sd + Sj;  ed = ed + nk;
        else;        He = He + Hj;  Se = Se + Sj;  ee = ee + nk;  end
    end

    if Sd + Se <= 0
        warning('momento_dobradica:vazio','nenhum painel de "%s" encontrado.', nome);
        return;
    end
    if Sd > 0
        H.cf_dir = Sd/ed;  H.Sf_dir = Sd;  H.H_dir = Hd;
        H.Ch_dir = Hd / (q * Sd * H.cf_dir);
    end
    if Se > 0
        H.cf_esq = Se/ee;  H.Sf_esq = Se;  H.H_esq = He;
        H.Ch_esq = He / (q * Se * H.cf_esq);
    end
    H.Ch_dif = H.Ch_dir - H.Ch_esq;
    cf_m = (Sd+Se)/max(ed+ee,eps);
    H.Ch_tot = (Hd+He) / (q * (Sd+Se) * cf_m);
    H.Ch_faixa = H.H_faixa ./ (q * H.Sf_faixa .* H.cf_faixa);

    CG = [obter(Aeronave,'Xcg',0), 0, obter(Aeronave,'Zcg',0)];
    Mt = [0 0 0];  Fz = 0;
    for i = 1:numel(P)
        A_no = P(i).P1 + 0.25*(P(i).P2 - P(i).P1);
        B_no = P(i).P4 + 0.25*(P(i).P3 - P(i).P4);
        Mt = Mt + cross((A_no+B_no)/2 - CG, P(i).Forca);
        Fz = Fz + P(i).Forca(3);
    end
    H.Cl_rol = Mt(1) / (q * Aeronave.Sref * Aeronave.Bref);
    H.Cn_guin= Mt(3) / (q * Aeronave.Sref * Aeronave.Bref);
    H.Cm_arf = Mt(2) / (q * Aeronave.Sref * Aeronave.Cref);
    H.CFz    = Fz    / (q * Aeronave.Sref);
end

function v = obter(s,nome,padrao)
    if isfield(s,nome) && ~isempty(s.(nome)); v = s.(nome); else; v = padrao; end
end
