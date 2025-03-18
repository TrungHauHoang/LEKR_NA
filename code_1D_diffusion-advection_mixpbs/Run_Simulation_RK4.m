function [Matrix_solution_RK4,vectork_rk4,vec_MVMs_4,vec_linear_combination_4,vec_linear_combination5vectors_4,vec_FLO_4,vec_err4_L2,vec_err4_max,vec_fetch4,vec_store4] = Run_Simulation_RK4(num_k_RK4,T,h,Jacobian,U0,U_ref,upper_bound_RK4,module_eig_max)

Matrix_solution_RK4 = zeros(length(U0),num_k_RK4);

vec_MVMs_4 = zeros(num_k_RK4,1);
vec_linear_combination_4 = zeros(num_k_RK4,1);
vec_linear_combination5vectors_4 = zeros(num_k_RK4,1);
vec_FLO_4 = zeros(num_k_RK4,1);
numsteps_arr4 = zeros(num_k_RK4,1);
vec_err4_L2 = zeros(num_k_RK4,1);
vec_err4_max = zeros(num_k_RK4,1);
vectork_rk4 = zeros(num_k_RK4,1);
vec_fetch4 = zeros(num_k_RK4,2);
vec_store4 = zeros(num_k_RK4,2);


for i = 1:num_k_RK4
    if i ~= 1
        numsteps_arr4(i) = numsteps_arr4(i-1)*2;
        step_size = T/numsteps_arr4(i);
    elseif i==1
        step_size = upper_bound_RK4/module_eig_max;
        numsteps_arr4(i) = ceil(T/step_size);
        step_size = T/numsteps_arr4(i);
    end
    for j = 1:50000
        t = 0:step_size:T;
        [ U1, MVMs1,linear_combination,Linear_combination5vectors, fetch, store, Operations ] = RK4_1D_count(Jacobian,t,U0);
        %% Calculate the error
        error_u = abs(U_ref-U1);
        vec_err4_L2(i) = (h*sum(error_u.^2))^0.5 
        vec_err4_max(i) = norm(U_ref-U1,'inf')
        if ( isnan(vec_err4_L2(i)) == 0 )
            vec_fetch4(i,1) = fetch(1);
            vec_store4(i,1) = store(1);
            vec_fetch4(i,2) = fetch(2);
            vec_store4(i,2) = store(2);
            vec_MVMs_4(i) = MVMs1;
            vec_linear_combination_4(i) = linear_combination;
            vec_linear_combination5vectors_4(i) = Linear_combination5vectors;
            vec_FLO_4(i) = Operations;
            Matrix_solution_RK4(:,i) = U1;
            break;
        else
            numsteps_arr4(i) = numsteps_arr4(i) + 1;
            step_size = T/numsteps_arr4(i);
        end
    end
    vectork_rk4(i) = step_size;
end

