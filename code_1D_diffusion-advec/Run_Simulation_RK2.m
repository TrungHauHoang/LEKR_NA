function [Matrix_solution_RK2,vectork_rk2,vec_MVMs_5,vec_linear_combination_5,vec_linear_combination3vectors_5,vec_err5_L2,vec_err5_max,vec_fetch5,vec_store5] = Run_Simulation_RK2(num_k_RK2,T,h,Jacobian,U0,U_ref,upper_bound_RK2,module_eig_max,dimension)


Matrix_solution_RK2 = zeros(length(U0),num_k_RK2);
vec_MVMs_5 = zeros(num_k_RK2,1);
vec_linear_combination_5 = zeros(num_k_RK2,1);
vec_linear_combination3vectors_5 = zeros(num_k_RK2,1);
numsteps_arr5 = zeros(num_k_RK2,1);
vec_err5_L2 = zeros(num_k_RK2,1);
vec_err5_max_L2 = zeros(num_k_RK2,1);
vec_err5_max = zeros(num_k_RK2,1);
vectork_rk2 = zeros(num_k_RK2,1);
vec_fetch5 = zeros(num_k_RK2,2);
vec_store5 = zeros(num_k_RK2,2);


for i = 1:num_k_RK2
    if i ~= 1
        numsteps_arr5(i) = numsteps_arr5(i-1)*2;
        step_size = T/numsteps_arr5(i);
    elseif i==1
        step_size = upper_bound_RK2/module_eig_max;
        numsteps_arr5(i) = ceil(T/step_size);
        step_size = T/numsteps_arr5(i);
    end
    
    for j = 1:50000
        t = 0:step_size:T;
        [ U1, MVMs1,linear_combination,Linear_combination3vectors, fetch, store  ] = RK2_1D_count(Jacobian,t,U0);
        error_u = abs(U_ref-U1);
        
            vec_err5_L2(i) = (h*sum(error_u.^2))^0.5
            vec_err5_max(i) = norm(U_ref-U1,'inf')
        
        if ( isnan(vec_err5_L2(i)) == 0  )
            vec_fetch5(i,1) = fetch(1);
            vec_store5(i,1) = store(1);
            vec_fetch5(i,2) = fetch(2);
            vec_store5(i,2) = store(2);
            vec_MVMs_5(i) = MVMs1;
            vec_linear_combination_5(i) = linear_combination;
            vec_linear_combination3vectors_5(i) = Linear_combination3vectors;
            Matrix_solution_RK2(:,i) = U1;
            break;
        else
            numsteps_arr5(i) = numsteps_arr5(i) + 1;
            step_size = T/numsteps_arr5(i);
        end
    end
    vectork_rk2(i) = step_size;
end

