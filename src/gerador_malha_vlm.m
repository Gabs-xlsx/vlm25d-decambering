function [Malha] = gerador_malha_vlm(Aeronave)
% GERADOR_MALHA_VLM  Reticulado do VLM: malha plana com normais de camber.

    Nx      = Aeronave.Malha.nx_corda;
    Ny_semi = Aeronave.Malha.ny_env;
    tipo_x  = obter_campo(Aeronave.Malha, 'espac_x', 1);
    tipo_y  = obter_campo(Aeronave.Malha, 'espac_y', 1);

    Sec_Dir = Aeronave.Asa.Secoes;
    Sec_Esq = Sec_Dir;
    for k = 1:numel(Sec_Esq); Sec_Esq(k).Y = -Sec_Dir(k).Y; end
    Secoes = [Sec_Esq(end:-1:2), Sec_Dir];

    Controles = obter_campo(Aeronave.Asa, 'Controles', []);

    quebras = [Secoes.Y];
    for c = 1:numel(Controles)
        quebras = [quebras, Controles(c).inicio_y, Controles(c).fim_y, ...
                           -Controles(c).inicio_y, -Controles(c).fim_y];
    end
    quebras = unique(round(quebras(:)' * 1e6) / 1e6);
    lim_y   = [min([Secoes.Y]), max([Secoes.Y])];
    quebras = quebras(quebras >= lim_y(1) - 1e-6 & quebras <= lim_y(2) + 1e-6);
    quebras = unique(min(max(quebras, lim_y(1)), lim_y(2)));
    if numel(quebras) < 2 || (quebras(end) - quebras(1)) <= 0
        error('gerador_malha_vlm:quebras', ...
              'quebras de envergadura degeneradas: %s (limites %s, %d secoes).', ...
              mat2str(quebras, 6), mat2str(lim_y, 6), numel(Secoes));
    end

    Ny_alvo = Ny_semi * 2;
    aperto  = obter_campo(Aeronave.Malha, 'aperto_ponta', 0.75);
    Y_nodes = espacamento_envergadura(quebras, Ny_alvo, tipo_y, aperto);
    Ny = numel(Y_nodes) - 1;

    nome_y = {'uniforme','global','por particao (LEGADO)'};
    if tipo_y == 1
        extra = sprintf(' %s, aperto de ponta %.2f', nome_y{2}, aperto);
    else
        extra = [' ' nome_y{min(max(tipo_y,0),2)+1}];
    end
    fprintf('    -> Malha VLM: %d x %d = %d paineis (espac x=%d, y=%d%s)\n', ...
            Nx, Ny, Nx*Ny, tipo_x, tipo_y, extra);

    x_hinge_ref = 0.75;
    if ~isempty(Controles)
        x_hinge_ref = mean(1.0 - [Controles.corda_c]);
    end
    Nx_front = min(max(round(Nx * x_hinge_ref), 1), Nx-1);
    Nx_flap  = Nx - Nx_front;

    S_frac = zeros(Nx+1, Ny+1);
    Par    = repmat(struct('xle',0,'zle',0,'c',1,'tor',0,'xcr',[],'zcr',[], ...
                           'xct',[],'zct',[],'f',0), 1, Ny+1);
    X_v  = zeros(Nx+1, Ny+1);  Y_v = zeros(Nx+1, Ny+1);
    Z_v  = zeros(Nx+1, Ny+1);
    Xc_v = zeros(Nx+1, Ny+1);  Zc_v = zeros(Nx+1, Ny+1);

    for j = 1:(Ny+1)
        y_loc = Y_nodes(j);

        x_hinge = x_hinge_ref;
        for c = 1:numel(Controles)
            if abs(y_loc) >= Controles(c).inicio_y - 1e-6 && ...
               abs(y_loc) <= Controles(c).fim_y    + 1e-6
                x_hinge = 1.0 - Controles(c).corda_c; break;
            end
        end
        tf = distribuir_espacamento(Nx_front, tipo_x);
        tl = distribuir_espacamento(Nx_flap,  tipo_x);
        X_frac = [x_hinge*tf(1:end-1), x_hinge + (1-x_hinge)*tl];

        k = find(y_loc >= [Secoes(1:end-1).Y] - 1e-6 & ...
                 y_loc <= [Secoes(2:end).Y]   + 1e-6, 1, 'first');
        if isempty(k); k = numel(Secoes) - 1; end
        f = (y_loc - Secoes(k).Y) / (Secoes(k+1).Y - Secoes(k).Y);
        f = min(max(f, 0), 1);

        x_le   = Secoes(k).X_LE  + f*(Secoes(k+1).X_LE  - Secoes(k).X_LE);
        z_le   = Secoes(k).Z_LE  + f*(Secoes(k+1).Z_LE  - Secoes(k).Z_LE);
        corda  = Secoes(k).Corda + f*(Secoes(k+1).Corda - Secoes(k).Corda);
        torcao = Secoes(k).Torcao+ f*(Secoes(k+1).Torcao- Secoes(k).Torcao);

        [xc_r, zc_r] = perfil_linha_media(Secoes(k).Perfil);
        [xc_t, zc_t] = perfil_linha_media(Secoes(k+1).Perfil);

        e_corda = [cosd(torcao), -sind(torcao)];
        e_perp  = [sind(torcao),  cosd(torcao)];

        Par(j) = struct('xle',x_le,'zle',z_le,'c',corda,'tor',torcao, ...
                        'xcr',xc_r,'zcr',zc_r,'xct',xc_t,'zct',zc_t,'f',f);
        S_frac(:,j) = X_frac(:);

        for i = 1:(Nx+1)
            s = X_frac(i);
            z_cam = interp1(xc_r, zc_r, s, 'linear', 'extrap') + ...
                f * (interp1(xc_t, zc_t, s, 'linear', 'extrap') - ...
                     interp1(xc_r, zc_r, s, 'linear', 'extrap'));

            X_v(i,j) = x_le + corda*s*e_corda(1);
            Z_v(i,j) = z_le + corda*s*e_corda(2);
            Y_v(i,j) = y_loc;

            Xc_v(i,j) = X_v(i,j) + corda*z_cam*e_perp(1);
            Zc_v(i,j) = Z_v(i,j) + corda*z_cam*e_perp(2);
        end
    end

    n_pan = Nx*Ny;
    Malha.Paineis(n_pan) = struct('ID',[],'i',[],'j',[],'P1',[],'P2',[],'P3',[], ...
        'P4',[],'VP',[],'CP',[],'n',[],'n0',[],'Area',[],'Corda',[], ...
        'is_control',[],'nome_controle',[]);
    idx = 1;
    for j = 1:Ny
        y_mid = abs((Y_v(1,j) + Y_v(1,j+1))/2);
        for i = 1:Nx
            P1 = [X_v(i,j),     Y_v(i,j),     Z_v(i,j)];
            P2 = [X_v(i+1,j),   Y_v(i+1,j),   Z_v(i+1,j)];
            P3 = [X_v(i+1,j+1), Y_v(i+1,j+1), Z_v(i+1,j+1)];
            P4 = [X_v(i,j+1),   Y_v(i,j+1),   Z_v(i,j+1)];

            P_ba = (P1+P4)/2;  P_bf = (P2+P3)/2;
            VP = P_ba + 0.25*(P_bf - P_ba);
            CP = P_ba + 0.75*(P_bf - P_ba);

            s_cp = S_frac(i,j) + 0.75*(S_frac(i+1,j) - S_frac(i,j));
            vx = tangente_camber(s_cp, j);
            Q1 = [Xc_v(i,j),     Y_v(i,j),     Zc_v(i,j)];
            Q2 = [Xc_v(i+1,j),   Y_v(i+1,j),   Zc_v(i+1,j)];
            Q3 = [Xc_v(i+1,j+1), Y_v(i+1,j+1), Zc_v(i+1,j+1)];
            Q4 = [Xc_v(i,j+1),   Y_v(i,j+1),   Zc_v(i,j+1)];
            vy = (Q4+Q3)/2 - (Q1+Q2)/2;
            nv = cross(vx, vy);  nv = nv / norm(nv);

            eh_ctrl = false;  nome_ctrl = '';
            if i > Nx_front
                for c = 1:numel(Controles)
                    if y_mid >= Controles(c).inicio_y - 1e-6 && ...
                       y_mid <= Controles(c).fim_y    + 1e-6
                        eh_ctrl = true; nome_ctrl = Controles(c).nome; break;
                    end
                end
            end

            Malha.Paineis(idx) = struct('ID',idx,'i',i,'j',j, ...
                'P1',P1,'P2',P2,'P3',P3,'P4',P4,'VP',VP,'CP',CP, ...
                'n',nv,'n0',nv, ...
                'Area',norm(cross(P3-P1, P4-P2))/2, ...
                'Corda',norm(P_bf - P_ba), ...
                'is_control',eh_ctrl,'nome_controle',nome_ctrl);
            idx = idx + 1;
        end
    end

    Malha.X_verts = X_v;  Malha.Y_verts = Y_v;  Malha.Z_verts = Z_v;
    Malha.Nx = Nx;  Malha.Ny = Ny;  Malha.Nx_front = Nx_front;
    Malha.Y_nodes = Y_nodes;

    function t = tangente_camber(s_cp, jj)
        P = Par(jj);
        h = 1e-5;
        sa = min(max(s_cp - h, 0), 1);
        sb = min(max(s_cp + h, 0), 1);
        [xa, za] = ponto_medio_local(P, sa);
        [xb, zb] = ponto_medio_local(P, sb);
        t = [xb - xa, 0, zb - za];
    end

    function [xx, zz] = ponto_medio_local(P, s)
        zr = interp1(P.xcr, P.zcr, s, 'linear', 'extrap');
        zt = interp1(P.xct, P.zct, s, 'linear', 'extrap');
        z_cam = zr + P.f*(zt - zr);
        ec = [cosd(P.tor), -sind(P.tor)];
        ep = [sind(P.tor),  cosd(P.tor)];
        xx = P.xle + P.c*s*ec(1) + P.c*z_cam*ep(1);
        zz = P.zle + P.c*s*ec(2) + P.c*z_cam*ep(2);
    end
end

function v = obter_campo(s, nome, padrao)
    if isfield(s, nome) && ~isempty(s.(nome)); v = s.(nome); else; v = padrao; end
end
