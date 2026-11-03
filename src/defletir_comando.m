function Malha = defletir_comando(Malha, nome, delta_deg, tipo)
% DEFLETIR_COMANDO  Gira fisicamente os paineis a jusante da dobradica.
%   Malha = defletir_comando(Malha, 'aileron', 15, 'assimetrico')

    if nargin < 4 || isempty(tipo); tipo = 'assimetrico'; end
    if abs(delta_deg) < 1e-12; return; end

    Nx = Malha.Nx;  Ny = Malha.Ny;
    if ~isfield(Malha,'Nx_front') || isempty(Malha.Nx_front)
        error('defletir_comando:sem_quebra', ...
              'Malha.Nx_front ausente -- a malha nao foi gerada com comando.');
    end
    nf = Malha.Nx_front;
    n_def = 0;

    for j = 1:Ny
        base = (j-1)*Nx;
        alvo = false;  y_ref = 0;
        for i = nf+1:Nx
            p = Malha.Paineis(base+i);
            if p.is_control && strcmpi(p.nome_controle, nome)
                alvo = true;  y_ref = p.CP(2);  break;
            end
        end
        if ~alvo; continue; end

        s = 1;
        if strcmpi(tipo,'assimetrico') && y_ref < 0; s = -1; end
        ang = s * deg2rad(delta_deg);

        pf = Malha.Paineis(base + nf + 1);
        A = pf.P1;  Bv = pf.P4;
        k = Bv - A;  nk = norm(k);
        if nk < 1e-12; continue; end
        k = k / nk;

        for i = nf+1:Nx
            p = Malha.Paineis(base+i);
            p.P1 = rodrigues(p.P1, A, k, ang);
            p.P2 = rodrigues(p.P2, A, k, ang);
            p.P3 = rodrigues(p.P3, A, k, ang);
            p.P4 = rodrigues(p.P4, A, k, ang);
            n0r  = p.n0*cos(ang) + cross(k,p.n0)*sin(ang) + k*(k*p.n0')*(1-cos(ang));
            p = atualizar_geometria_painel(p);
            p.n0 = n0r / norm(n0r);
            p.n  = p.n0;
            Malha.Paineis(base+i) = p;
            n_def = n_def + 1;
        end
    end

    Malha.deflexao.(matlab.lang.makeValidName(nome)) = delta_deg;
    Malha.deflexao_tipo = tipo;
    if n_def == 0
        warning('defletir_comando:nada', ...
                'nenhum painel marcado como "%s" -- nada foi defletido.', nome);
    end
end

function q = rodrigues(p, A, k, ang)
    v = p - A;
    q = A + v*cos(ang) + cross(k,v)*sin(ang) + k*(k*v')*(1-cos(ang));
end
