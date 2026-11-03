function painel = atualizar_geometria_painel(painel)
% =========================================================================

    P1 = painel.P1;
    P2 = painel.P2;
    P3 = painel.P3;
    P4 = painel.P4;

    P_front = (P1 + P4) / 2;
    P_back  = (P2 + P3) / 2;

    painel.VP = P_front + 0.25 * (P_back - P_front);
    painel.CP = P_front + 0.75 * (P_back - P_front);

    vec_x = P_back - P_front;
    vec_y = (P4+P3)/2 - (P1+P2)/2;
    n_novo = cross(vec_x, vec_y);
    painel.n = n_novo / norm(n_novo);

    painel.Area = norm(cross(P3-P1, P4-P2)) / 2;

end
