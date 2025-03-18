
%%========================================================================

function [Jacobian,MVMs] = Jacobian_matrix(minus_A_diff_x,A_diffusion,kappa)
Jacobian = minus_A_diff_x + kappa*A_diffusion;
MVMs = 1;
end

