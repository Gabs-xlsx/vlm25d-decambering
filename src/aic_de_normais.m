function A = aic_de_normais(C, Paineis)
% AIC_DE_NORMAIS  Contrai o tensor de influencia com as normais atuais.

    N  = numel(Paineis);
    nn = reshape([Paineis.n], 3, N)';
    A  = C(:,:,1) .* nn(:,1) + C(:,:,2) .* nn(:,2) + C(:,:,3) .* nn(:,3);
end
