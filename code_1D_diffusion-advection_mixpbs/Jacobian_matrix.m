
%%========================================================================

function [Jacobian,MVMs,fetch,store,Operations] = Jacobian_matrix(minus_A_diff_x,A_diffusion,kappa)

fetch = zeros(2,1);
store = zeros(2,1);

Jacobian = minus_A_diff_x + kappa.*A_diffusion;

MVMs = 1;
fetch(1) = 0;
store(1) = 0;
fetch(2) = 2;
store(2) = 1;
Operations = 2*nnz(A_diffusion);
end

