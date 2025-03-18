function [Matrix_solution_Leja,U1,vec_MVMs_3,vec_inner_product_3,vec_Scalar_multiplication_3,vec_linear_combination_3,vec_err3_L2,vec_err3_max,vec_fetch3,vec_store3] = Run_Simulation_Leja(dof,T,h,numsteps,index_work,k,kappa,U0,U0_ori,U_ref,CFL,type_eigenvalue,num_k,Jacobian,tolerence)


vec_MVMs_3 = zeros(num_k,1);
vec_inner_product_3 = zeros(num_k,1);
vec_Scalar_multiplication_3 = zeros(num_k,1);
vec_linear_combination_3 = zeros(num_k,1);
vec_err3_L2 = zeros(num_k,1);
vec_err3_max = zeros(num_k,1);
vec_fetch3 = zeros(num_k,2);
vec_store3 = zeros(num_k,2);
CFL_ori = CFL;

k_ori = k;
if (strcmp(type_eigenvalue,'REAL') == 1)
    for i = index_work:num_k-index_work+1
            if (i>=8) 
                i = i+2;
            end
            
        spectrum_end = min(real(eig(full(Jacobian))));
        eigenvalue_min = k*spectrum_end;
        eigenvalue_max = 0;
            
        for j =1:numsteps
            norm_U0 = norm(U0);
            [U1,MVMs3,inner_product,Scalar_multiplication,linear_combination,fetch3,store3,Operations3] = real_leja_degree(k,Jacobian,U0/norm_U0,eigenvalue_min,eigenvalue_max,tolerence(i),0);
            %Calculate the costs
            U1 = norm_U0*U1;
            U0 = U1;
            vec_MVMs_3(i) = vec_MVMs_3(i)+ MVMs3;
            vec_inner_product_3(i) = vec_inner_product_3(i) + inner_product + 1;
            vec_Scalar_multiplication_3(i) = vec_Scalar_multiplication_3(i) + Scalar_multiplication + 2;
            vec_linear_combination_3(i) = vec_linear_combination_3(i) + linear_combination;
            vec_fetch3(i,1) = vec_fetch3(i,1) + fetch3(1) + 1;
            vec_store3(i,1) = vec_store3(i,1) + store3(1) + 1;
            
        end
        %% Calculate the error
        error_u = abs(U_ref-U1);
        vec_err3_L2(i) =  (h*sum(error_u.^2))^0.5
        vec_err3_max(i) = norm(U_ref-U1,'inf')
        
        %% Returning to original, to compute the next time step
        
        U0 = U0_ori;
        k = k/2;
        numsteps = ceil(T/k);
        k = T/numsteps;
        
        if i==7
            U0 = U0_ori;
            k = k_ori;
            numsteps = ceil(T/k);
            k = T/numsteps;
        end
        
    end
elseif (strcmp(type_eigenvalue,'COMPLEX') == 1)
    Matrix_solution_Leja = {};
    
    for i = index_work:num_k

        spectrum_adv = max(abs(imag(eig(full(Jacobian)))));
        eigenvalue_min = -spectrum_adv;
        eigenvalue_max = spectrum_adv;
        Matrix_solution = [];
        for j =1:numsteps
            norm_U0 = norm(U0);
            [U1,MVMs3,inner_product,Scalar_multiplication,linear_combination,fetch3,store3] = cplx_leja_degree(k,Jacobian,U0/norm_U0,eigenvalue_min,eigenvalue_max,tolerence(i),0,0);
            U1 = norm_U0*U1;
            U0 = U1;
            Matrix_solution = [Matrix_solution;U0];
           
            vec_MVMs_3(i) = vec_MVMs_3(i)+ MVMs3;
            vec_inner_product_3(i) = vec_inner_product_3(i) + inner_product + 1;
            vec_Scalar_multiplication_3(i) = vec_Scalar_multiplication_3(i) + Scalar_multiplication + 2;
            vec_linear_combination_3(i) = vec_linear_combination_3(i) + linear_combination;
            vec_fetch3(i,1) = vec_fetch3(i,1) + fetch3(1) + 1;
            vec_store3(i,1) = vec_store3(i,1) + store3(1) + 1;
            
        end
        Matrix_solution_Leja{i} = Matrix_solution;
        %% Calculate the error
        error_u = abs(U_ref-U1);
        vec_err3_L2(i) =  (h*sum(error_u.^2))^0.5
        vec_err3_max(i) = norm(U_ref-U1,'inf')
        
        %% Returning to original, to compute the next time step
        U0 = U0_ori;
        k = k/2;
        numsteps = ceil(T/k);
        k = T/numsteps;
        if i==5
            U0 = U0_ori;
            k = k_ori;
            numsteps = ceil(T/k);
            k = T/numsteps;
        end
    end
end