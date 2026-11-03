function C = montar_tensor_influencia(Malha, x_infinito)
% MONTAR_TENSOR_INFLUENCIA  Tensor de influencia aerodinamica, N x N x 3.

    if nargin < 2 || isempty(x_infinito); x_infinito = 1e6; end

    Paineis = Malha.Paineis;
    N = numel(Paineis);
    CP = reshape([Paineis.CP], 3, N)';

    C = zeros(N, N, 3);
    for j = 1:N
        P1 = Paineis(j).P1;  P2 = Paineis(j).P2;
        P3 = Paineis(j).P3;  P4 = Paineis(j).P4;

        A_j = P1 + 0.25*(P2 - P1);
        B_j = P4 + 0.25*(P3 - P4);
        A_inf = [x_infinito, A_j(2), A_j(3)];
        B_inf = [x_infinito, B_j(2), B_j(3)];

        V =   biot_savart_multi(A_inf, A_j, CP) ...
            + biot_savart_multi(A_j,   B_j, CP) ...
            + biot_savart_multi(B_j, B_inf, CP);

        C(:, j, :) = reshape(V, N, 1, 3);
    end
end
