function [Matrix_solution_Krylov,vectork,evaluateF_1,vec_MVMs_1,vec_inner_product_1,vec_Scalar_multiplication_1,vec_linear_combination_1,vec_err1_L2,vec_err1_max_L2,vec_err1_max,vec_fetch1,vec_store1] = Run_Simulation_Krylov(T,theta,h,numsteps,index_work,k,internalPoints,nu,cases,size_rho,U0,U0_ori,U_ref,CFL,num_k)

Matrix_solution_Krylov = zeros(length(U0),num_k);
vec_MVMs_1 = zeros(num_k,1);
vec_inner_product_1 = zeros(num_k,1);
vec_linear_combination_1 = zeros(num_k,1);
vec_Scalar_multiplication_1 = zeros(num_k,1);
evaluateF_1 = zeros(num_k,1);
vec_err1_L2 = zeros(num_k,1);
vec_err1_max_L2 = zeros(num_k,1);
vec_err1_max = zeros(num_k,1);
vectork = zeros(num_k,1);
vec_fetch1 = zeros(num_k,2);
vec_store1 = zeros(num_k,2);
dof = internalPoints+1;

logspace_vector_kry = logspace(log10( 1e-3 ), log10( 0.5*(1e-7+1e-6)), num_k);

for i = index_work:num_k
    
    for j =1:numsteps

        F_Un = vectorF_2D(dof,cases,nu,size_rho,U0);
        normF = norm(F_Un);
        [w1,MVMs3,inner_product,Scalar_multiplication,linear_combination,fetch3,store3] = Krylov_sparse(dof,cases,nu,size_rho,k,U0,F_Un/normF,110,1,logspace_vector_kry(i));
        U1 = U0 + normF*k*w1;
        U0 = U1;
        
        evaluateF_1(i) = evaluateF_1(i) + 1;
        vec_MVMs_1(i) = vec_MVMs_1(i) + MVMs3;
        vec_inner_product_1(i) = vec_inner_product_1(i) + inner_product + 1;
        vec_Scalar_multiplication_1(i) = vec_Scalar_multiplication_1(i) + Scalar_multiplication + 1;
        vec_linear_combination_1(i) = vec_linear_combination_1(i) + linear_combination + 1;
        vec_fetch1(i,1) = vec_fetch1(i,1)  + fetch3(1) + 3;
        vec_store1(i,1) = vec_store1(i,1)  + store3(1) + 3;
    end
    Matrix_solution_Krylov(:,i) = U1;
    %% Calculate the error
    error_rho = abs(U_ref(1:size_rho)-U1(1:size_rho));
    error_u = abs(U_ref(size_rho+1:2*size_rho)-U1(size_rho+1:2*size_rho));
    error_v = abs(U_ref(2*size_rho+1:end)-U1(2*size_rho+1:end));
    vec_err1_L2(i) = h*norm(error_rho,2) + h*norm(error_u,2) + h*norm(error_v,2)
    vec_err1_max_L2(i) = max([h*norm(error_rho,2),h*norm(error_u,2),h*norm(error_v,2)])
    vec_err1_max(i) = norm(U_ref-U1,'inf')
    
    %% Returning to original, to compute the next time step
    U0 = U0_ori;
    vectork(i) = k;
    CFL = CFL/2;
    k = CFL*min([0.5*h,(h^2/(8*nu)),(h/(theta))]);
    numsteps = ceil(T/k);
    k = T/numsteps;
    
end

for i=index_work-1:-1:1
    vectork(i) = vectork(i+1)*2;
end

