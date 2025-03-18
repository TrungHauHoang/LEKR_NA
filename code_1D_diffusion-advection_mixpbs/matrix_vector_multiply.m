function v = matrix_vector_multiply(kappa,h,u)
    n = length(u);
    v = zeros(n, 1);

    % Apply finite difference stencil for second-order derivative
    for i = 1:n
        if i == 1
            v(i) = 2*u(i) - u(i+1);  % Homogeneous Dirichlet boundary condition
        elseif i == n
            v(i) = 2*u(i) - u(i-1);  % Homogeneous Dirichlet boundary condition
        else
            v(i) = -u(i-1) + 2*u(i) - u(i+1);
        end
    end
end