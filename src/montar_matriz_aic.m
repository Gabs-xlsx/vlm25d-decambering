function [A, C] = montar_matriz_aic(Malha)
% MONTAR_MATRIZ_AIC  Compatibilidade: monta o tensor e ja devolve o AIC.

    C = montar_tensor_influencia(Malha);
    A = aic_de_normais(C, Malha.Paineis);
end
