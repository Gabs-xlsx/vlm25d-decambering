%EXEMPLO  Asa retangular AR 6 com aileron, perfil S1223, Re 2,5e5 (~2 min na 1a execucao).
%   Calcula a solucao potencial e a acoplada 2,5D para alpha = 5 graus e
%   aileron defletido em +10 graus, e imprime C_L e C_h da superficie direita.
%   Requer xfoil.exe nesta pasta ou no PATH do sistema.

addpath(fullfile(fileparts(mfilename('fullpath')), 'src'));

Re = 250e3;  xtr = [1 1];  ncrit = 9;          % transicao livre
alpha = 5;  delta = 10;                         % graus
V = 3.77;  rho = 1.225;

% ---- geometria: semi-envergadura 3 m, corda 1 m, aileron em 57% a 100% da semi-envergadura
AIL = struct('nome','aileron','tipo','assimetrico','inicio_eta',0.60/1.05,'fim_eta',1.00,'corda_c',0.25);
G = struct('nome','AR6','semi_env',3.0,'corda_raiz',1.0,'corda_ponta',1.0, ...
           'perfil','s1223.dat','superficies',AIL);
Ae = montar_aeronave(G);
Ae.Malha.nx_corda = 32;  Ae.Malha.ny_env = 20;
Ae.CasosVoo(1).V_inf = V;  Ae.CasosVoo(1).rho = rho;

% ---- polares 2D (XFOIL + Viterna) e familia defletida
Pe = extrapolar_viterna(polar_estoque('s1223.dat', Re, -12:0.5:20, struct('xtr',xtr,'ncrit',ncrit)));
Pol = repmat(Pe, 1, numel(Ae.Asa.Secoes));
xh = 1 - AIL.corda_c;
FAM = familia_estoque('s1223.dat', Re, [-delta delta], xh, numel(Ae.Asa.Secoes), xtr, ncrit);
Fam = struct('delta',{delta,-delta},'Polares',{[],[]});
Fam(1).Polares = FAM(abs([FAM.delta] - delta) < 1e-9).Polares;
Fam(2).Polares = FAM(abs([FAM.delta] + delta) < 1e-9).Polares;

% ---- malha e deflexao
M0 = gerador_malha_vlm(Ae);
Ms = defletir_comando(M0, 'aileron', delta, AIL.tipo);
Ms = polares_por_faixa(Ms, Pol, Ae, Fam, struct('nome',''));

% ---- solucao potencial
Cs = Ae.CasosVoo(1);  Cs.alpha = alpha;
Rp = resolver_vlm(Ms, montar_matriz_aic(Ms), Cs, Ae);
Hp = momento_dobradica(Ms, Rp, Ae, 'aileron', V, rho);

% ---- solucao acoplada 2,5D (configuracao adotada no trabalho)
opts = struct('esquema',2,'x2',0.80,'usar_tara',true,'eta_lib',0.95, ...
              'lim_estado',deg2rad(40),'modo_raiz','max','jac_consistente',true, ...
              'damp',0.5,'tol',1e-3,'max_iter',25,'verbose',false, ...
              'modo_liberacao','segmento','n_extrap',2,'n_borda',2);
[Ra, MD] = resolver_vlm_25D(Ms, Cs, Ae, Pol, opts);
Ha = momento_dobradica(MD, Ra, Ae, 'aileron', V, rho);

fprintf('\nalpha = %g graus, delta = %+g graus\n', alpha, delta);
fprintf('  potencial : C_L = %.4f   C_h = %+.4f\n', Rp.CL, Hp.Ch_dir);
fprintf('  acoplado  : C_L = %.4f   C_h = %+.4f   (convergiu: %d)\n', Ra.CL, Ha.Ch_dir, Ra.convergiu);
