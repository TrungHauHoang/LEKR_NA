function [largest_eigen_value] = Power_iteration(dof,cases,nu,size_rho,u)
    % Parameters
    % ----------
    % u                       : Input state variable(s)
    % RHS_function            : RHS function
    %
    % Returns
    % -------
    % largest_eigen_value     : Largest eigenvalue (within 2% accuracy)
    % rhs_calls               : Number of RHS calls (3*ii)
    
    tol = 1e-4;               % 2% tolerance
    niters = 1000;            % Max. number of iterations
    eigenvalue_ii_1 = 0;      % Eigenvalue at ii-1
    vector = ones(size(u));   % Initial estimate of eigenvector
    
    for ii = 1:niters
        
        % Compute new eigenvector
        eigenvector = Jacobian_actvec_approx(dof,cases,nu,size_rho,u, vector);
        
        % Norm of eigenvector = eigenvalue
        eigenvalue = norm(eigenvector);
        
        % Normalize eigenvector to eigenvalue; new estimate of eigenvector
        vector = eigenvector / eigenvalue;
        
        % Check convergence for eigenvalues (eigenvalues converge faster than eigenvectors)
        if abs(eigenvalue - eigenvalue_ii_1) <= tol * eigenvalue
            largest_eigen_value = eigenvalue;
            break;
        end
        
        % This value becomes the previous one
        eigenvalue_ii_1 = eigenvalue;
        
    end

end

