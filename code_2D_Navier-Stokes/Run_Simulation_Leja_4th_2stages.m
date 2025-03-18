function [Matrix_solution_Leja,evaluateF_3,U1,vec_MVMs_3,vec_inner_product_3,vec_Scalar_multiplication_3,vec_linear_combination_3,vec_err3_L2,vec_err3_max_L2,vec_err3_max,vec_fetch3,vec_store3] = Run_Simulation_Leja_4th_2stages(T,theta,h,numsteps,index_work,k,internalPoints,nu,cases,size_rho,U0,U0_ori,U_ref,CFL,num_k)

Matrix_solution_Leja = zeros(length(U0),num_k);
vec_MVMs_3 = zeros(num_k,1);
vec_inner_product_3 = zeros(num_k,1);
vec_Scalar_multiplication_3 = zeros(num_k,1);
vec_linear_combination_3 = zeros(num_k,1);
evaluateF_3 = zeros(num_k,1);
vec_err3_L2 = zeros(num_k,1);
vec_err3_max_L2 = zeros(num_k,1);
vec_err3_max = zeros(num_k,1);
vec_fetch3 = zeros(num_k,2);
vec_store3 = zeros(num_k,2);
dof = internalPoints+1;

logspace_vector_Leja = logspace(log10(1e-9), log10(1e-18), num_k);
logspace_vector_Leja1 = logspace(log10(1e2), log10(1e10), num_k);

for i = index_work:num_k
    
    for j =1:numsteps
        %Calculate the Rosenbrock-Euler methods
        F_Un = vectorF_2D(dof,cases,nu,size_rho,U0);
        if (j == 1) || (mod(j,50) == 0)
            
            adv_max = Power_iteration(dof,cases,nu,size_rho,U0);
            adv_max = 1.1*adv_max
            adv_min = -adv_max
            
        end
        g_n = F_Un - Jacobian_matrix_actvec(dof,cases,nu,size_rho,U0,U0);
        %% For computation of Un2
        
        normF = norm(F_Un);
        [w1,MVMs3_Un2,inner_product1,Scalar_multiplication_Un2,linear_combination_Un2,fetch3_Un2,store3_Un2] = cplx_leja_sparse_degree(dof,cases,nu,size_rho,(3/4)*k,U0,adv_min,adv_max,F_Un/normF,400,1,0,logspace_vector_Leja(i)*logspace_vector_Leja1(i));
        Un2 = U0 + (3/4)*normF*k*w1;
        
        %% For computation of Un+1
        F_Un2 = vectorF_2D(dof,cases,nu,size_rho,Un2);
        g_Un2 = F_Un2 - Jacobian_matrix_actvec(dof,cases,nu,size_rho,U0,Un2);
        V0_4th = [U0; zeros(3,1)];
        V0_4th(end) = 1;
        norm_V0_4th = norm(V0_4th);
        [w1,MVMs3_Un,inner_product2,Scalar_multiplication_Un,linear_combination_Un,fetch3_Un,store3_Un] = cplx_leja_sparse_degree(dof,cases,nu,size_rho,k,U0,adv_min,adv_max,V0_4th/norm_V0_4th,400,0,1,logspace_vector_Leja(i),((1/k)^2)*(32/9)*(g_Un2-g_n),zeros(size(g_n)),g_n);
        w1 = norm_V0_4th*w1;
        U1 = w1(1:end-3);
        evaluateF_3(i) = evaluateF_3(i) + 2;
        vec_MVMs_3(i) = vec_MVMs_3(i) + MVMs3_Un2 + MVMs3_Un + 32;
        vec_inner_product_3(i) = vec_inner_product_3(i)+ inner_product1 + inner_product2 + 2;
        vec_Scalar_multiplication_3(i) = vec_Scalar_multiplication_3(i) + Scalar_multiplication_Un2 + Scalar_multiplication_Un + 4;
        vec_linear_combination_3(i) = vec_linear_combination_3(i) + linear_combination_Un2 + linear_combination_Un + 1;
        vec_fetch3(i,1) = vec_fetch3(i,1) + fetch3_Un2(1) + fetch3_Un(1) + 27;
        vec_store3(i,1) = vec_store3(i,1) + store3_Un2(1) + store3_Un(1) + 30;
        U0 = U1;
        
    end
    Matrix_solution_Leja(:,i) = U1;
    %% Calculate the error
    error_rho = abs(U_ref(1:size_rho)-U1(1:size_rho));
    error_u = abs(U_ref(size_rho+1:2*size_rho)-U1(size_rho+1:2*size_rho));
    error_v = abs(U_ref(2*size_rho+1:end)-U1(2*size_rho+1:end));
    vec_err3_L2(i) = h*(sum(error_rho.^2))^0.5 + h*(sum(error_u.^2))^0.5 + h*(sum(error_v.^2))^0.5
    vec_err3_max_L2(i) = max([h*(sum(error_rho.^2))^0.5,h*(sum(error_u.^2))^0.5 ,h*(sum(error_v.^2))^0.5])
    vec_err3_max(i) = norm(U_ref-U1,'inf')
    
    %% Returning to original, to compute the next time step
    U0 = U0_ori;
    CFL = CFL/2;
    k = CFL*min([0.5*h,(h^2/(8*nu)),(h/(theta))]);
    numsteps = ceil(T/k);
    k = T/numsteps;
    
end


