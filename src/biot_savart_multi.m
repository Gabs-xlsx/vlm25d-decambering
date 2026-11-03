function V = biot_savart_multi(P1, P2, PTS)
% BIOT_SAVART_MULTI  Versao vetorizada de biot_savart_segmento.

    r1 = PTS - P1;
    r2 = PTS - P2;
    r1xr2 = cross(r1, r2, 2);
    den = sum(r1xr2.^2, 2);

    n1 = sqrt(sum(r1.^2, 2));
    n2 = sqrt(sum(r2.^2, 2));

    num = sum((r1 - r2) .* (r1./n1 - r2./n2), 2);

    K = num ./ (4*pi*den);
    K(den < 1e-10 | n1 < 1e-10 | n2 < 1e-10) = 0;

    V = K .* r1xr2;
end
