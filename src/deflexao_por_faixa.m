function [dj, info] = deflexao_por_faixa(Malha, nome)
%DEFLEXAO_POR_FAIXA  Le a deflexao LOCAL de cada faixa direto da geometria.
%   [dj, info] = deflexao_por_faixa(Malha, 'aileron')

    Nx = Malha.Nx;  Ny = Malha.Ny;
    dj = zeros(Ny,1);
    if ~isfield(Malha,'Nx_front') || isempty(Malha.Nx_front)
        info = struct('carimbo',NaN,'discorda',NaN,'n_faixas',0);
        return;
    end
    nf = Malha.Nx_front;
    if nf < 1 || nf >= Nx
        info = struct('carimbo',NaN,'discorda',NaN,'n_faixas',0);
        return;
    end

    n_faixas = 0;
    for j = 1:Ny
        base = (j-1)*Nx;
        alvo = false;
        for i = nf+1:Nx
            p = Malha.Paineis(base+i);
            if p.is_control && (isempty(nome) || strcmpi(p.nome_controle, nome))
                alvo = true; break;
            end
        end
        if ~alvo; continue; end
        n_faixas = n_faixas + 1;

        pa = Malha.Paineis(base + nf);
        pb = Malha.Paineis(base + nf + 1);

        k = pb.P4 - pb.P1;  nk = norm(k);
        if nk < 1e-12; continue; end
        k = k / nk;

        u = 0.5*((pa.P2 - pa.P1) + (pa.P3 - pa.P4));
        v = 0.5*((pb.P2 - pb.P1) + (pb.P3 - pb.P4));
        nu = norm(u);  nv = norm(v);
        if nu < 1e-12 || nv < 1e-12; continue; end
        u = u/nu;  v = v/nv;

        kp = k - dot(k,u)*u;  nkp = norm(kp);
        if nkp < 1e-12; continue; end
        kp = kp / nkp;
        dj(j) = rad2deg(atan2(dot(cross(u,v), kp), dot(u,v)));
    end

    dj(abs(dj) < 1e-9) = 0;

    car = NaN;
    if isfield(Malha,'deflexao') && ~isempty(nome)
        cn = matlab.lang.makeValidName(nome);
        if isfield(Malha.deflexao, cn); car = Malha.deflexao.(cn); end
    end
    disc = NaN;
    if isfinite(car)
        sel = false(Ny,1);
        for j = 1:Ny
            base = (j-1)*Nx;
            for i = nf+1:Nx
                p = Malha.Paineis(base+i);
                if p.is_control && (isempty(nome) || strcmpi(p.nome_controle,nome))
                    sel(j) = true; break;
                end
            end
        end
        if any(sel); disc = max(abs(abs(dj(sel)) - abs(car))); end
        if isfinite(disc) && disc > 0.05
            warning('deflexao_por_faixa:discorda', ...
                ['a geometria e o carimbo Malha.deflexao discordam em ate %.3f deg ' ...
                 '(carimbo %.2f). Vale o MEDIDO.'], disc, car);
        end
    end

    info = struct('carimbo',car,'discorda',disc,'n_faixas',n_faixas, ...
                  'deltas_unicos',unique(round(dj*1e4)/1e4));
end
