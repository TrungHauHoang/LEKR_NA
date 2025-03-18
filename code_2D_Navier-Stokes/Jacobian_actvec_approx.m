function Jacobian_vector = Jacobian_actvec_approx(dof,cases,nu,size_rho,u, v)
    % Jacobian_vector = (RHS(u + epsilon*v) - RHS(u))/epsilon
    
    % epsilon is normalized to norm(u)
    epsilon = 1e-7 * norm(u);

    % J(u) * v = (RHS(u + epsilon*v) - RHS(u))/epsilon
    Jacobian_vector = (vectorF_2D(dof,cases,nu,size_rho,u + (epsilon * v)) - vectorF_2D(dof,cases,nu,size_rho,u)) / epsilon;
end


