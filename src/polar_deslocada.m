function Pf = polar_deslocada(P, delta_deg, x_dob)
%POLAR_DESLOCADA  Polar do perfil COM flape, por teoria de perfil fino.
%   Pf = polar_deslocada(P, 10, 0.80)

    if nargin < 3 || isempty(x_dob); x_dob = 0.80; end
    T = teoria_flape(x_dob);
    d_a  = T.dalphaL0(delta_deg);
    d_cm = T.m2 * tan(deg2rad(delta_deg));

    Pf = P;
    Pf.alpha = P.alpha(:) + d_a;
    Pf.Cm    = P.Cm(:)    + d_cm;
    if isfield(P,'alpha_pedido')  && ~isempty(P.alpha_pedido);  Pf.alpha_pedido  = P.alpha_pedido(:)  + d_a; end
    if isfield(P,'alpha_faltando')&& ~isempty(P.alpha_faltando); Pf.alpha_faltando = P.alpha_faltando(:) + d_a; end
    if isfield(P,'alpha_min_real') && ~isempty(P.alpha_min_real); Pf.alpha_min_real = P.alpha_min_real + d_a; end
    if isfield(P,'alpha_max_real') && ~isempty(P.alpha_max_real); Pf.alpha_max_real = P.alpha_max_real + d_a; end

    Pf.Origem = 'teoria de perfil fino (deslocamento)';
    if isfield(P,'Nome')
        Pf.Nome = sprintf('%s | flape %+.1f deg (teoria)', P.Nome, delta_deg);
    end
    Pf.flape = struct('delta_deg',delta_deg,'x_dob',x_dob, ...
                      'd_alphaL0_deg',d_a,'d_Cm',d_cm,'fonte','teoria');
end
