
clear all
close all
format long

%% Parameters to change
ax = 0.0;
bx = 1.0;
dimension = 2;
cases = 2; 
case_initial = 8; 
internalPoints = 159; 
h = (bx-ax)/(internalPoints+1);

theta = 1/12;
C = 1e-6*((1/theta)*(1/h));

% theta = 1/6;
% C = 1e-4*((1/theta)*(1/h));

T = 1/theta;
CFL = 320;
alpha = 1;
%% Compute initial data
nu = theta*h*C;
CFL_ori = CFL;
[U0,size_rho,dof,index_work_array] = initial2D(ax,bx,internalPoints,h,theta,cases,case_initial,C,CFL_ori);
U0_ori = U0;
sizeU0 = length(U0);
k = CFL*min([0.5*h,(h^2/(8*nu)),(h/(theta))]);
numsteps_ori = ceil(T/k);
k_ori = T/numsteps_ori;
CFL_actual = (CFL_ori*k_ori)/k;
num_k_Kry = 8;
num_k_Che = 8;
num_k_Leja = 8;

%% Create Jacobian matrix

Jacobian = Jacobian_matrix(dof,cases,nu,size_rho,U0_ori);
Jacobian_ori = Jacobian;
eig_max = eigs(Jacobian_ori,1,'largestabs','Tolerance',1e-5);
module_eig_max = abs(eig_max);
if isreal(eig_max) == 1
    upper_bound_RK4 = 2.7853;
    upper_bound_RK2 = 2;
    type_eigenvalue = 'REAL';
else
    upper_bound_RK4 = 2.7;
    upper_bound_RK2 = 0.40;%0.44;
    type_eigenvalue = 'COMPLEX';
end

%% Reference solution

[U_ref] = index_2D(C,CFL,internalPoints,case_initial,theta,h);

%% =========== 1/ Using Krylov subspace methods ============================

CFL = CFL_ori;
index_work = index_work_array(1);
CFL = CFL/(2^(index_work-1));
k = CFL*min([0.5*h,(h^2/(8*nu)),(h/(theta))]);
numsteps = ceil(T/k);
k = T/numsteps;
tStart = tic;
[Matrix_solution_Krylov,vectork,evaluateF_1,vec_MVMs_1,vec_inner_product_1,vec_Scalar_multiplication_1,vec_linear_combination_1,vec_err1_L2,vec_err1_max_L2,vec_err1_max,vec_fetch1,vec_store1] = Run_Simulation_Krylov(T,theta,h,numsteps,index_work,k,internalPoints,nu,cases,size_rho,U0,U0_ori,U_ref,CFL,num_k_Kry);
tEnd = toc(tStart)

n_Jacobianw1 = (vec_MVMs_1)/16;
vec_MVMs_1 = vec_MVMs_1 + 10*evaluateF_1;
Total_Fetch_Krylov = 18*n_Jacobianw1 + 12*evaluateF_1 + 6*vec_linear_combination_1 + 3*vec_Scalar_multiplication_1 + 6*alpha*vec_inner_product_1 + vec_fetch1(:,1) ;
Total_Store_Krylov = 3*n_Jacobianw1 + 12*evaluateF_1 + 3*vec_linear_combination_1 + 3*vec_Scalar_multiplication_1 + vec_store1(:,1);
Total_MOPS_Krylov = 21*n_Jacobianw1 + 12*evaluateF_1 + 9*vec_linear_combination_1 + 6*vec_Scalar_multiplication_1 + 6*alpha*vec_inner_product_1 + vec_fetch1(:,1)  + vec_store1(:,1);

Total_MOPS_Krylov_1 = 21*n_Jacobianw1 + 12*evaluateF_1 + 9*vec_linear_combination_1 + 6*vec_Scalar_multiplication_1 + 6*alpha*10*vec_inner_product_1 + vec_fetch1(:,1)  + vec_store1(:,1);

%% =========== 2/ Using Leja methods ===============================================

CFL = CFL_ori; 
index_work = index_work_array(2);
CFL = CFL/(2^(index_work-1));
k = CFL*min([0.5*h,(h^2/(8*nu)),(h/(theta))]);
numsteps = ceil(T/k);
k = T/numsteps;
tStart = tic;
[Matrix_solution_Leja,evaluateF_2,U_Leja_T,vec_MVMs_2,vec_inner_product_2,vec_Scalar_multiplication_2,vec_linear_combination_2,vec_err2_L2,vec_err2_max_L2,vec_err2_max,vec_fetch2,vec_store2] = Run_Simulation_Leja(T,theta,h,numsteps,index_work,k,internalPoints,nu,cases,size_rho,U0,U0_ori,U_ref,CFL,num_k_Leja);
tEnd = toc(tStart)

n_Jacobianw3 = (vec_MVMs_2)/16;
vec_MVMs_2 = vec_MVMs_2 + 10*evaluateF_2;
Total_Fetch_Leja = 18*n_Jacobianw3 + 12*evaluateF_2 + 6*vec_linear_combination_2 + 3*vec_Scalar_multiplication_2 + 6*alpha*vec_inner_product_2 + vec_fetch2(:,1);
Total_Store_Leja = 3*n_Jacobianw3 + 12*evaluateF_2 + 3*vec_linear_combination_2 + 3*vec_Scalar_multiplication_2 + + vec_store2(:,1) ;
Total_MOPS_Leja = 21*n_Jacobianw3 + 12*evaluateF_2 + 9*vec_linear_combination_2 + 6*vec_Scalar_multiplication_2 + 9*alpha*vec_inner_product_2  + vec_fetch2(:,1) + vec_store2(:,1);

%% 4th order Leja two stages
CFL = CFL_ori; 
index_work = index_work_array(3);
CFL = CFL/(2^(index_work-1));
k = CFL*min([0.5*h,(h^2/(8*nu)),(h/(theta))]);
numsteps = ceil(T/k);
k = T/numsteps;
tStart = tic;
[Matrix_solution_Leja_4th2S,evaluateF_3_4th2S,U_Leja_T_4th2S,vec_MVMs_3_4th2S,vec_inner_product_3_4th2S,vec_Scalar_multiplication_3_4th2S,vec_linear_combination_3_4th2S,vec_err3_L2_4th2S,vec_err3_max_L2_4th2S,vec_err3_max_4th2S,vec_fetch3_4th2S,vec_store3_4th2S] = Run_Simulation_Leja_4th_2stages(T,theta,h,numsteps,index_work,k,internalPoints,nu,cases,size_rho,U0,U0_ori,U_ref,CFL,num_k_Leja);
tEnd = toc(tStart)

n_Jacobianw3_4th2S = (vec_MVMs_3_4th2S)/16;
vec_MVMs_3_4th2S = vec_MVMs_3_4th2S + 10*evaluateF_3_4th2S;
Total_Fetch_Leja_4th2S = 18*n_Jacobianw3_4th2S + 12*evaluateF_3_4th2S + 6*vec_linear_combination_3_4th2S + 3*vec_Scalar_multiplication_3_4th2S + 6*alpha*vec_inner_product_3_4th2S + vec_fetch3_4th2S(:,1);
Total_Store_Leja_4th2S = 3*n_Jacobianw3_4th2S + 12*evaluateF_3_4th2S + 3*vec_linear_combination_3_4th2S + 3*vec_Scalar_multiplication_3_4th2S + + vec_store3_4th2S(:,1) ;
Total_MOPS_Leja_4th2S = 21*n_Jacobianw3_4th2S + 12*evaluateF_3_4th2S + 9*vec_linear_combination_3_4th2S + 6*vec_Scalar_multiplication_3_4th2S + 9*alpha*vec_inner_product_3_4th2S + vec_fetch3_4th2S(:,1) + vec_store3_4th2S(:,1);

%% RK4

U0 = U0_ori;
num_k_RK4 = 5;
step_sizerk4 = logspace(log10(upper_bound_RK4/module_eig_max), log10(upper_bound_RK4/(module_eig_max*16)), num_k_RK4);
numsteps_arr4 = ceil(T./step_sizerk4);
step_sizerk4 = T./numsteps_arr4;

tStart = tic;
[Matrix_solution_RK4,vectork_rk4,evaluateF_4,vec_MVMs_4,vec_inner_product_4,vec_linear_combination_4,vec_linear_combination5vectors_4,vec_err4_L2,vec_err4_max_L2,vec_err4_max,vec_fetch4,vec_store4] = Run_Simulation_RK4_cost(num_k_RK4,T,h,internalPoints,nu,cases,size_rho,U0,U_ref,step_sizerk4);
tEnd = toc(tStart)

vec_MVMs_4 = vec_MVMs_4 + 10*evaluateF_4;
Total_Fetch_RK4 = 12*evaluateF_4 + 6*vec_linear_combination_4 + 15*vec_linear_combination5vectors_4  + vec_fetch4(:,1);
Total_Store_RK4 = 12*evaluateF_4 + 3*vec_linear_combination_4 + 3*vec_linear_combination5vectors_4  + vec_store4(:,1);
Total_MOPS_RK4 = 12*evaluateF_4 + 9*vec_linear_combination_4 + 18*vec_linear_combination5vectors_4  + vec_fetch4(:,1) + vec_store4(:,1);

%% RK2

U0 = U0_ori;
num_k_RK2 = 5;
step_sizerk2 = logspace(log10(upper_bound_RK2/module_eig_max), log10(upper_bound_RK2/(module_eig_max*16)), num_k_RK2);
numsteps_arr2 = ceil(T./step_sizerk2);
step_sizerk2 = T./numsteps_arr2;

tStart = tic;
[Matrix_solution_RK2,vectork_rk2,evaluateF_5,vec_MVMs_5,vec_inner_product_5,vec_linear_combination_5,vec_linear_combination3vectors_5,vec_err5_L2,vec_err5_max_L2,vec_err5_max,vec_fetch5,vec_store5] = Run_Simulation_RK2_cost(num_k_RK2,T,h,internalPoints,nu,cases,size_rho,U0,U_ref,step_sizerk2);
tEnd = toc(tStart)

vec_MVMs_5 = vec_MVMs_5 + 10*evaluateF_5;
Total_Fetch_RK2 = 12*evaluateF_5 + 6*vec_linear_combination_5 + 9*vec_linear_combination3vectors_5  + vec_fetch5(:,1);
Total_Store_RK2 = 12*evaluateF_5 + 3*vec_linear_combination_5 + 3*vec_linear_combination3vectors_5 + vec_store5(:,1);
Total_MOPS_RK2 = 12*evaluateF_5 + 9*vec_linear_combination_5 + 12*vec_linear_combination3vectors_5  + vec_fetch5(:,1) + vec_store5(:,1);

%% Back to original

k = k_ori;
numsteps = numsteps_ori;
CFL = CFL_ori;
U0 = U0_ori;

%% =========== compare ====================================================
%% Memory operations
h_plot = figure;

loglog( Total_MOPS_Krylov, vec_err1_L2,'-g^','LineWidth',2,'MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',5);
hold on
loglog( Total_MOPS_Krylov_1, vec_err1_L2,'--g^','LineWidth',2,'MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',5);
loglog(Total_MOPS_Leja,vec_err2_L2,'-mo','LineWidth',2,'MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',5);
loglog(Total_MOPS_Leja_4th2S,vec_err3_L2_4th2S,'--mo','LineWidth',2,'MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',5);
loglog(Total_MOPS_RK4, vec_err4_L2,'-s','Color',[0.6350, 0.0780, 0.1840],'LineWidth',2,'MarkerEdgeColor',[0.6350, 0.0780, 0.1840],'MarkerFaceColor',[0.6350, 0.0780, 0.1840],'MarkerSize',5);
loglog(Total_MOPS_RK2, vec_err5_L2,'-bd','LineWidth',2,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',5);
ylim([-inf 10]);
hold off
set(h_plot,'Units','Inches');
pos = get(h_plot,'Position');
set(h_plot,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h_plot,'filename','-dpdf','-r0')
set(gca,'FontSize',15)
set(gca,'TickLength',[0.03,0])
xlabel('Total number of memory operations','FontSize', 17);
ylabel('Errors in discrete L^2 norm','FontSize', 17);
legend({'Krylov \zeta=1','Krylov \zeta=10','Leja','Leja4th2S','RK4','RK2'},'location','best','NumColumns',2,'FontSize',13);

