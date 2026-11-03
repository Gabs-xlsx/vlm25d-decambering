function t = distribuir_espacamento(n, tipo)
% DISTRIBUIR_ESPACAMENTO  Vetor de n+1 nos em [0,1].

    if nargin < 2 || isempty(tipo); tipo = 1; end
    if tipo == 0
        t = linspace(0, 1, n+1)';
    else
        t = 0.5 * (1 - cos(linspace(0, pi, n+1)'));
    end
    t = t(:)';
end
