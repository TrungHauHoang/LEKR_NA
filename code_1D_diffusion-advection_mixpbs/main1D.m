
clear all
close all
format short
%% Parameters to change

ax = 0.0;
bx = 1.0;
dimension = 1;
cases_boundary = 1; 
case_initial = 1; 
internalPoints = 159;
h = (bx-ax)/(internalPoints+1);
T = 1;
C = 2;
CFL = 1280;
alpha = 1;

%% Compute initial data
CFL_ori = CFL;
[U0,dof,index_work_array] = initial_1D(ax,bx,internalPoints,cases_boundary,case_initial,CFL);
U0_ori = U0;
sizeU0 = length(U0);
k_ori = 0.1;
num_k = 14;

%% Create Jacobian matrix

mesh_x = linspace(ax, bx, 159); % Generate 1000 points between 0 and 1
kappa = 33/5120 + (31/5120)*(tanh(20*(mesh_x - 0.8)));

A_nabla = lap1d_nabla(dof,cases_boundary);
minus_A_diff_x = -A_nabla;
A_diffusion = lap1d(dof,cases_boundary);
Jacobian = Jacobian_matrix(minus_A_diff_x,A_diffusion,kappa);

%% Reference solution
tStart = tic;
eig_max = eigs(Jacobian,1,'largestabs','Tolerance',1e-5)
module_eig_max = abs(eig_max);

if isreal(eig_max) == 1
    upper_bound_RK4 = 2.6;
    upper_bound_RK2 = 1.9;
    type_eigenvalue = 'REAL';
else
    upper_bound_RK4 = 2.7;
    upper_bound_RK2 = 1.05;
    type_eigenvalue = 'COMPLEX';
end

step_temp = (upper_bound_RK4/module_eig_max)/64;
N_tilda = ceil(T/step_temp);
step = T/N_tilda;
t = 0:step:T;
[ U_ref, MVMs_Uref_RK4, fetch_Uref_RK4, store_Uref_RK4, Operations_Uref_RK4 ] = RK4_1D_count(Jacobian,t,U0);

tEnd = toc(tStart)
%% =========== 1/ Using Krylov subspace methods ============================

k = k_ori;
index_work = index_work_array(1);
numsteps = ceil(T/k);
tStart = tic;
[Matrix_solution_Krylov,vectork,vec_MVMs_1,vec_inner_product_1,vec_Scalar_multiplication_1,vec_linear_combination_1,vec_FLO_1,vec_err1_L2,vec_err1_max,vec_fetch1,vec_store1] = Run_Simulation_Krylov(T,h,numsteps,index_work,k,U0,U0_ori,U_ref,num_k,type_eigenvalue,Jacobian);
tEnd = toc(tStart)
Total_Fetch_Krylov = vec_MVMs_1 + 2*vec_linear_combination_1 + vec_Scalar_multiplication_1 + 2*alpha*vec_inner_product_1 + vec_fetch1(:,1);
Total_Store_Krylov = vec_MVMs_1 + vec_linear_combination_1 + vec_Scalar_multiplication_1 + vec_store1(:,1);
Total_MOPS_Krylov = 2*vec_MVMs_1 + 3*vec_linear_combination_1 + 2*vec_Scalar_multiplication_1 + 2*alpha*vec_inner_product_1 + vec_fetch1(:,1) + vec_store1(:,1);

%% =========== 2/ Using Leja methods ===============================================

k = k_ori;
index_work = index_work_array(2);
numsteps = ceil(T/k);
tStart = tic;
[Matrix_solution_Leja,U_Leja_T,vec_MVMs_2,vec_inner_product_2,vec_Scalar_multiplication_2,vec_linear_combination_2,vec_FLO_2,vec_err2_L2,vec_err2_max,vec_fetch2,vec_store2] =  Run_Simulation_Leja(T,h,sizeU0,numsteps,index_work,k,U0,U0_ori,U_ref,type_eigenvalue,num_k,Jacobian);
tEnd = toc(tStart)

Total_Fetch_Leja = vec_MVMs_2 + 2*vec_linear_combination_2 + vec_Scalar_multiplication_2 + 2*alpha*vec_inner_product_2 + vec_fetch2(:,1);
Total_Store_Leja = vec_MVMs_2 + vec_linear_combination_2 + vec_Scalar_multiplication_2 + vec_store2(:,1);
Total_MOPS_Leja = 2*vec_MVMs_2 + 3*vec_linear_combination_2 + 2*vec_Scalar_multiplication_2 + 2*alpha*vec_inner_product_2 + vec_fetch2(:,1) + vec_store2(:,1);

%% RK4

num_k_RK4 = 5;
tStart = tic;
[Matrix_solution_RK4,vectork_rk4,vec_MVMs_3,vec_linear_combination_3,vec_linear_combination5vectors_3,vec_FLO_3,vec_err3_L2,vec_err3_max,vec_fetch3,vec_store3] = Run_Simulation_RK4(num_k_RK4,T,h,Jacobian,U0,U_ref,upper_bound_RK4,module_eig_max);
tEnd = toc(tStart)
Total_Fetch_RK4 = vec_MVMs_3 + 2*vec_linear_combination_3 + 5*vec_linear_combination5vectors_3  + vec_fetch3(:,1);
Total_Store_RK4 = vec_MVMs_3 + vec_linear_combination_3 + vec_linear_combination5vectors_3  + vec_store3(:,1);
Total_MOPS_RK4 = 2*vec_MVMs_3 + 3*vec_linear_combination_3 + 6*vec_linear_combination5vectors_3  + vec_fetch3(:,1) + vec_store3(:,1) ;

%% RK2
num_k_RK2 = 5;
tStart = tic;
[Matrix_solution_RK2,vectork_rk2,vec_MVMs_4,vec_linear_combination_4,vec_linear_combination3vectors_4,vec_FLO_4,vec_err4_L2,vec_err4_max,vec_fetch4,vec_store4] = Run_Simulation_RK2(num_k_RK2,T,h,Jacobian,U0,U_ref,upper_bound_RK2,module_eig_max);
tEnd = toc(tStart)

Total_Fetch_RK2 = vec_MVMs_4 + 2*vec_linear_combination_4 + 3*vec_linear_combination3vectors_4  + vec_fetch4(:,1);
Total_Store_RK2 = vec_MVMs_4 + vec_linear_combination_4 + vec_linear_combination3vectors_4 + vec_store4(:,1);
Total_MOPS_RK2 = 2*vec_MVMs_4 + 3*vec_linear_combination_4 + 4*vec_linear_combination3vectors_4  + vec_fetch4(:,1) + vec_store4(:,1);

%% Back to original

k = k_ori;
CFL = CFL_ori;
U0 = U0_ori;
num_k1 = num_k/2;

%% Error in discrete L2 norm
h_plot = figure;

loglog( vectork(index_work_array(1):num_k1-1), vec_err1_L2(index_work_array(1):num_k1-1),'-.g^','LineWidth',2,'MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',5);
hold on
loglog(vectork(index_work_array(2)+1:num_k1-1),vec_err2_L2(index_work_array(2)+1:num_k1-1),'-.mo','LineWidth',2,'MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',5);
loglog(vectork_rk4, vec_err3_L2,'-s','Color',[0.6350, 0.0780, 0.1840],'LineWidth',2,'MarkerEdgeColor',[0.6350, 0.0780, 0.1840],'MarkerFaceColor',[0.6350, 0.0780, 0.1840],'MarkerSize',5);
loglog(vectork_rk2, vec_err4_L2,'-bd','LineWidth',2,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',5);
loglog( vectork(num_k1+index_work_array(1):end-1), vec_err1_L2(num_k1+index_work_array(1):end-1),'-g^','LineWidth',2,'MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',5);
loglog(vectork(num_k1+index_work_array(2):end-1), vec_err2_L2(num_k1+index_work_array(2):end-1),'-mo','LineWidth',2,'MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',5);
loglog(vectork(1:num_k1), vectork(1:num_k1).^2,'-k.','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',0.001);
loglog(vectork(1:num_k1), vectork(1:num_k1).^4,'--k.','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',0.001);
ylim([-inf 10]);
xlim([1e-4/2 2]);
hold off
set(h_plot,'Units','Inches');
pos = get(h_plot,'Position');
set(h_plot,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h_plot,'filename','-dpdf','-r0')

set(gca,'FontSize',13)
set(gca,'TickLength',[0.03,0])
xlabel('Time step size \tau','FontSize', 17);
ylabel('Errors in discrete L^2 norm','FontSize', 17);
legend({'Krylov 10^{-4}','Leja 10^{-4} ','RK4','RK2','Krylov 10^{-7}','Leja 10^{-7}','f(\tau) = c\tau^2','g(\tau) = c\tau^4'},'location','best','NumColumns',2,'FontSize', 13);

%% 
h_plot = figure;
 
loglog( vectork(index_work_array(1):num_k1-1), Total_MOPS_Krylov(index_work_array(1):num_k1-1),'-.g^','LineWidth',2,'MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',5);
hold on
loglog( vectork(index_work_array(2)+1:num_k1-1), Total_MOPS_Leja(index_work_array(2)+1:num_k1-1),'-.mo','LineWidth',2,'MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',5);
loglog(vectork_rk4, Total_MOPS_RK4(:,1),'-s','Color',[0.6350, 0.0780, 0.1840],'LineWidth',2,'MarkerEdgeColor',[0.6350, 0.0780, 0.1840],'MarkerFaceColor',[0.6350, 0.0780, 0.1840],'MarkerSize',5);
loglog( vectork(num_k1+index_work_array(1):end-1), Total_MOPS_Krylov(num_k1+index_work_array(1):end-1),'-g^','LineWidth',2,'MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',5);
loglog( vectork(num_k1+index_work_array(2)+1:end-1), Total_MOPS_Leja(num_k1+index_work_array(2)+1:end-1),'-mo','LineWidth',2,'MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',5);
loglog(vectork_rk2, Total_MOPS_RK2(:,1),'-bd','LineWidth',2,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',5);

hold off
set(h_plot,'Units','Inches');
pos = get(h_plot,'Position');
set(h_plot,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h_plot,'filename','-dpdf','-r0')
xlim([1e-4/2 2]);
set(gca,'FontSize',13)
set(gca,'TickLength',[0.03,0])
xlabel('Time step size \tau','FontSize', 17);
ylabel('Total number of memory operations','FontSize', 17);
legend({'Krylov 10^{-4}','Leja 10^{-4}','RK4','Krylov 10^{-7}','Leja 10^{-7}','RK2'},'location','best','NumColumns',2,'FontSize', 13);

%% spectrum of Jacobian matrix
full_eig  = eig((kappa.*A_diffusion + A_nabla));
h_plot = figure;
plot(real(full_eig),imag(full_eig),'*')
set(h_plot,'Units','Inches');
pos = get(h_plot,'Position');
set(h_plot,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h_plot,'filename','-dpdf','-r0')
xlabel('Re(z)');
ylabel('Im(z)');
set(gca,'FontSize',15)
set(gca,'TickLength',[0.025,0])
xlim([-1400 0])
ylim([-200 200])
legend('mixed problem')





