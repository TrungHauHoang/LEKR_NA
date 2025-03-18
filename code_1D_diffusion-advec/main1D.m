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
alpha = 1;

%% advection-dominated problem
C = 1/8;
CFL = 40;
%% diffusion-dominated problem
% C = 2;
% CFL = 1280;

%% Compute initial data
kappa = h*C;
CFL_ori = CFL;
[U0,dof,index_work_array,num_k,tolerence_kry,tolerence_leja] = initial_1D(ax,bx,internalPoints,cases_boundary,case_initial,CFL);
U0_ori = U0;
sizeU0 = length(U0);
k = CFL*min([h,(h^2/(4*kappa))]);
numsteps_ori = ceil(T/k);
k_ori = T/numsteps_ori;
CFL_actual = (CFL_ori*k_ori)/k;

%% Create Jacobian matrix

A_nabla = lap1d_nabla(dof,cases_boundary);
minus_A_diff_x = -A_nabla;
A_diffusion = lap1d(dof,cases_boundary);
Jacobian = Jacobian_matrix(minus_A_diff_x,A_diffusion,kappa);
Jacobian_ori = Jacobian;

%% Reference solution
tStart = tic;
eig_max = eigs(Jacobian_ori,1,'largestabs','Tolerance',1e-5)
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

CFL = CFL_ori;
index_work = index_work_array(1);
CFL = CFL/(2^(index_work-1));
k = CFL*min([h,(h^2/(4*kappa))]);
numsteps = ceil(T/k);
k = T/numsteps
tStart = tic;
[Matrix_solution_Krylov,vectork,vec_MVMs_1,vec_inner_product_1,vec_Scalar_multiplication_1,vec_linear_combination_1,vec_err1_L2,vec_err1_max,vec_fetch1,vec_store1] = Run_Simulation_Krylov(dof,T,h,numsteps,index_work,k,kappa,U0,U0_ori,U_ref,CFL,num_k,type_eigenvalue,Jacobian,tolerence_kry);
tEnd = toc(tStart)
Total_Fetch_Krylov = vec_MVMs_1 + 2*vec_linear_combination_1 + vec_Scalar_multiplication_1 + 2*alpha*vec_inner_product_1 + vec_fetch1(:,1) ;
Total_Store_Krylov = vec_MVMs_1 + vec_linear_combination_1 + vec_Scalar_multiplication_1 + vec_store1(:,1);
Total_MOPS_Krylov = 2*vec_MVMs_1 + 3*vec_linear_combination_1 + 2*vec_Scalar_multiplication_1 + 2*alpha*vec_inner_product_1 + vec_fetch1(:,1) + vec_store1(:,1);

%% =========== 2/ Using Leja methods ===============================================

CFL = CFL_ori; 
index_work = index_work_array(2);
CFL = CFL/(2^(index_work-1));
k = CFL*min([h,(h^2/(4*kappa))]);
numsteps = ceil(T/k);
k = T/numsteps;
tStart = tic;
[Matrix_solution_Leja,U_Leja_T,vec_MVMs_2,vec_inner_product_2,vec_Scalar_multiplication_2,vec_linear_combination_2,vec_err2_L2,vec_err2_max,vec_fetch2,vec_store2] =  Run_Simulation_Leja(dof,T,h,numsteps,index_work,k,kappa,U0,U0_ori,U_ref,CFL,type_eigenvalue,num_k,Jacobian,tolerence_leja);
tEnd = toc(tStart)

Total_Fetch_Leja = vec_MVMs_2 + 2*vec_linear_combination_2 + vec_Scalar_multiplication_2 + 2*alpha*vec_inner_product_2 + vec_fetch2(:,1);
Total_Store_Leja = vec_MVMs_2 + vec_linear_combination_2 + vec_Scalar_multiplication_2 + vec_store2(:,1);
Total_MOPS_Leja = 2*vec_MVMs_2 + 3*vec_linear_combination_2 + 2*vec_Scalar_multiplication_2 + 2*alpha*vec_inner_product_2 + vec_fetch2(:,1) + vec_store2(:,1);

%% RK4

num_k_RK4 = 5;
tStart = tic;
[Matrix_solution_RK4,vectork_rk4,vec_MVMs_3,vec_linear_combination_3,vec_linear_combination5vectors_3,vec_err3_L2,vec_err3_max,vec_fetch3,vec_store3] = Run_Simulation_RK4(num_k_RK4,T,h,Jacobian,U0,U_ref,upper_bound_RK4,module_eig_max,dimension);
tEnd = toc(tStart)
Total_Fetch_RK4 = vec_MVMs_3 + 2*vec_linear_combination_3 + 5*vec_linear_combination5vectors_3  + vec_fetch3(:,1);
Total_Store_RK4 = vec_MVMs_3 + vec_linear_combination_3 + vec_linear_combination5vectors_3  + vec_store3(:,1);
Total_MOPS_RK4 = 2*vec_MVMs_3 + 3*vec_linear_combination_3 + 6*vec_linear_combination5vectors_3  + vec_fetch3(:,1) + vec_store3(:,1);

%% RK2
num_k_RK2 = 5;
tStart = tic;
[Matrix_solution_RK2,vectork_rk2,vec_MVMs_4,vec_linear_combination_4,vec_linear_combination3vectors_4,vec_err4_L2,vec_err4_max,vec_fetch4,vec_store4] = Run_Simulation_RK2(num_k_RK2,T,h,Jacobian,U0,U_ref,upper_bound_RK2,module_eig_max,dimension);
tEnd = toc(tStart)
Total_Fetch_RK2 = vec_MVMs_4 + 2*vec_linear_combination_4 + 3*vec_linear_combination3vectors_4  + vec_fetch4(:,1);
Total_Store_RK2 = vec_MVMs_4 + vec_linear_combination_4 + vec_linear_combination3vectors_4 + vec_store4(:,1);
Total_MOPS_RK2 = 2*vec_MVMs_4 + 3*vec_linear_combination_4 + 4*vec_linear_combination3vectors_4  + vec_fetch4(:,1) + vec_store4(:,1);

%% Back to original

k = k_ori;
numsteps = numsteps_ori;
CFL = CFL_ori;
U0 = U0_ori;
num_k1 = num_k/2;


%% Error in discrete L2 norm
h_plot = figure;

loglog( vectork(index_work_array(1):num_k1), vec_err1_L2(index_work_array(1):num_k1),'-.g^','LineWidth',2,'MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',5);
hold on
loglog(vectork(index_work_array(2):num_k1),vec_err2_L2(index_work_array(2):num_k1),'-.mo','LineWidth',2,'MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',5);
loglog(vectork_rk4, vec_err3_L2,'-s','Color',[0.6350, 0.0780, 0.1840],'LineWidth',2,'MarkerEdgeColor',[0.6350, 0.0780, 0.1840],'MarkerFaceColor',[0.6350, 0.0780, 0.1840],'MarkerSize',5);
loglog(vectork_rk2, vec_err4_L2,'-bd','LineWidth',2,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',5);
loglog( vectork(num_k1+index_work_array(1):end), vec_err1_L2(num_k1+index_work_array(1):end),'-g^','LineWidth',2,'MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',5);
loglog(vectork(num_k1+index_work_array(2):end), vec_err2_L2(num_k1+index_work_array(2):end),'-mo','LineWidth',2,'MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',5);
loglog(vectork(1:num_k1), vectork(1:num_k1).^2,'-k.','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',0.001);
loglog(vectork(1:num_k1), vectork(1:num_k1).^4,'--k.','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',0.001);
ylim([-inf 10]);
xlim([1e-4/1.5 2]);
hold off
set(h_plot,'Units','Inches');
pos = get(h_plot,'Position');
set(h_plot,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h_plot,'filename','-dpdf','-r0')

set(gca,'FontSize',13.1)
set(gca,'TickLength',[0.03,0])
xlabel('Time step size \tau','FontSize', 17);
ylabel('Errors in discrete L^2 norm','FontSize', 17);
legend({'Krylov 10^{-4}','Leja 10^{-4}','RK4','RK2','Krylov 10^{-7}','Leja 10^{-7}','f(\tau) = c\tau^2','g(\tau) = c\tau^4'},'location','best','NumColumns',2,'FontSize', 13);

%% 
h_plot = figure;
 
loglog( vectork(index_work_array(1):num_k1), Total_MOPS_Krylov(index_work_array(1):num_k1),'-.g^','LineWidth',2,'MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',5);
hold on
loglog( vectork(index_work_array(2):num_k1), Total_MOPS_Leja(index_work_array(2):num_k1),'-.mo','LineWidth',2,'MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',5);
loglog(vectork_rk4, Total_MOPS_RK4(:,1),'-s','Color',[0.6350, 0.0780, 0.1840],'LineWidth',2,'MarkerEdgeColor',[0.6350, 0.0780, 0.1840],'MarkerFaceColor',[0.6350, 0.0780, 0.1840],'MarkerSize',5);
loglog( vectork(num_k1+index_work_array(1):end), Total_MOPS_Krylov(num_k1+index_work_array(1):end),'-g^','LineWidth',2,'MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',5);
loglog( vectork(num_k1+index_work_array(2):end), Total_MOPS_Leja(num_k1+index_work_array(2):end),'-mo','LineWidth',2,'MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',5);
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
legend({'Krylov 10^{-4}','Leja 10^{-4}','RK4','Krylov 10^{-7}','Leja 10^{-7}','RK2'},'location','best','NumColumns',1,'FontSize', 13);


%% spectrum of Jacobian matrix
matrix_temp = (1/320)*A_diffusion - A_nabla;
matrix_temp2 = 0.0125*A_diffusion - A_nabla;
full_eig  = eig( matrix_temp);
full_eig2  = eig( matrix_temp2);
h_plot = figure;
plot(real(full_eig2),imag(full_eig2),'*')
hold on
plot(real(full_eig),imag(full_eig),'*')
set(h_plot,'Units','Inches');
pos = get(h_plot,'Position');
set(h_plot,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h_plot,'filename','-dpdf','-r0')
xlabel('Re(z)');
ylabel('Im(z)');
set(gca,'FontSize',15)
set(gca,'TickLength',[0.025,0])
xlim([-1400 400])
ylim([-200 200])
legend('diffusion-dominated','advection-dominated')

%% spectrum of Jacobian matrix
x = ax+h:h:bx-h;
step_temp = (upper_bound_RK4/module_eig_max)/32;
N_tilda = ceil(0.25/step_temp);
step = 0.25/N_tilda;
t = 0:step:0.25;
[ U_ref1, MVMs_Uref_RK4, fetch_Uref_RK4, store_Uref_RK4, Operations_Uref_RK4 ] = RK4_1D_count(Jacobian,t,U0);

step_temp = (upper_bound_RK4/module_eig_max)/32;
N_tilda = ceil(0.5/step_temp);
step = 0.5/N_tilda;
t = 0:step:0.5;
[ U_ref2, MVMs_Uref_RK4, fetch_Uref_RK4, store_Uref_RK4, Operations_Uref_RK4 ] = RK4_1D_count(Jacobian,t,U0);

figure
plot(x, U_ref1);
figure
plot(x, U_ref2);
figure
plot(x, U_ref);

% U = Matrix_solution_Leja{1};
% for i = 1:numsteps*2
%     figure
%     plot(x, U((i-1)*159+1:i*159));
% end





