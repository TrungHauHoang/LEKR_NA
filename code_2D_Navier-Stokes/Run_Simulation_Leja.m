function [Matrix_solution_Leja,evaluateF_3,U1,vec_MVMs_3,vec_inner_product_3,vec_Scalar_multiplication_3,vec_linear_combination_3,vec_err3_L2,vec_err3_max_L2,vec_err3_max,vec_fetch3,vec_store3] = Run_Simulation_Leja(T,theta,h,numsteps,index_work,k,internalPoints,nu,cases,size_rho,U0,U0_ori,U_ref,CFL,num_k)

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
logspace_vector_Leja = logspace(log10(1e-7), log10(1e-10), num_k);

for i = index_work:num_k
    for j =1:numsteps
        F_Un = vectorF_2D(dof,cases,nu,size_rho,U0);
        if (j == 1) || (mod(j,50) == 0)
            
            adv_max = Power_iteration(dof,cases,nu,size_rho,U0);
            adv_max = 1.1*adv_max
            adv_min = -adv_max
            
        end
        normF = norm(F_Un);
        [w1,MVMs3,inner_product,Scalar_multiplication,linear_combination,fetch3,store3] = cplx_leja_sparse_degree(dof,cases,nu,size_rho,k,U0,adv_min,adv_max,F_Un/normF,400,1,0,logspace_vector_Leja(i));
        U1 = U0 + normF*k*w1;
        
        evaluateF_3(i) = evaluateF_3(i) + 1;
        vec_MVMs_3(i) = vec_MVMs_3(i) + MVMs3;
        vec_inner_product_3(i) = vec_inner_product_3(i) + 1 + inner_product;
        vec_Scalar_multiplication_3(i) = vec_Scalar_multiplication_3(i) + Scalar_multiplication + 1;
        vec_linear_combination_3(i) = vec_linear_combination_3(i) + linear_combination + 1;
        vec_fetch3(i,1) = vec_fetch3(i,1)  + fetch3(1) + 3;
        vec_store3(i,1) = vec_store3(i,1)  + store3(1) + 3;
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

