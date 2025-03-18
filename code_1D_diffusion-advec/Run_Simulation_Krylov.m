function [Matrix_solution_Krylov,vectork,vec_MVMs_1,vec_inner_product_1,vec_Scalar_multiplication_1,vec_linear_combination_1,vec_err1_L2,vec_err1_max,vec_fetch1,vec_store1] = Run_Simulation_Krylov(dof,T,h,numsteps,index_work,k,kappa,U0,U0_ori,U_ref,CFL,num_k,type_eigenvalue,Jacobian,tolerence)


Matrix_solution_Krylov = zeros(length(U0),num_k);
vec_MVMs_1 = zeros(num_k,1);
vec_inner_product_1 = zeros(num_k,1);
vec_Scalar_multiplication_1 = zeros(num_k,1);
vec_linear_combination_1 = zeros(num_k,1);
vec_err1_L2 = zeros(num_k,1);
vec_err1_max = zeros(num_k,1);
vectork = zeros(num_k,1);
vec_fetch1 = zeros(num_k);
vec_store1 = zeros(num_k);

n = size(Jacobian,1);
k_ori = k;

for i = index_work:num_k
    scaled_Jacobian = k*Jacobian;
    for j =1:numsteps
        norm_U0 = norm(U0);
        [U1,MVMs3,inner_product,Scalar_multiplication,linear_combination,fetch3,store3] = Krylov_sparse(scaled_Jacobian,U0/norm_U0,300,n,0,tolerence(i));
        U1 = norm_U0*U1;
        U0 = U1;
        %Calculate the costs
        vec_MVMs_1(i) = vec_MVMs_1(i) + MVMs3;
        vec_inner_product_1(i) = vec_inner_product_1(i) + inner_product + 1;
        vec_Scalar_multiplication_1(i) = vec_Scalar_multiplication_1(i) + Scalar_multiplication + 2;
        vec_linear_combination_1(i) = vec_linear_combination_1(i) + linear_combination;
        vec_fetch1(i) = vec_fetch1(i) + fetch3(1) + 1;
        vec_store1(i) = vec_store1(i) + store3(1) + 1;
    end
    %% Calculate the error
    error_u = abs(U_ref-U1);
    vec_err1_L2(i) = (h*sum(error_u.^2))^0.5
    vec_err1_max(i) = norm(U_ref-U1,'inf')
    
    %% Returning to original, to compute the next time step
    
    U0 = U0_ori;
    vectork(i) = k;
    k = k/2;
    numsteps = ceil(T/k);
    k = T/numsteps;
    
    if (i==7 && strcmp(type_eigenvalue,'REAL') == 1) ||  (i==5 && strcmp(type_eigenvalue,'COMPLEX') == 1)
        U0 = U0_ori;
        k = k_ori;
        numsteps = ceil(T/k);
        k = T/numsteps;
    end
end

end
