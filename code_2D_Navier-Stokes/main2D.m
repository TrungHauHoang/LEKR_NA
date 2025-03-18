
clear all
close all
format short


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
num_k = 8;

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
    upper_bound_RK2 = 0.4;
    type_eigenvalue = 'COMPLEX';
end

%% Reference solution
tStart = tic;

eig_max = eigs(Jacobian_ori,1,'largestabs','Tolerance',1e-5)

% step_temp = (upper_bound_RK4/module_eig_max)/64;
% N_tilda = ceil(T/step_temp);
% step = T/N_tilda;
% t = 0:step:T;
% [ U_ref, MVMs_Uref_RK4, fetch_Uref_RK4, store_Uref_RK4, Operations_Uref_RK4 ] = RK4_2D_count(internalPoints,nu,cases,size_rho,t,U0);

[U_ref] = index_2D(C,CFL,internalPoints,case_initial,theta,h);

tend = toc(tStart)
%% =========== 1/ Using Krylov subspace methods ============================

CFL = CFL_ori;
index_work = index_work_array(1);
CFL = CFL/(2^(index_work-1));
k = CFL*min([0.5*h,(h^2/(8*nu)),(h/(theta))]);
numsteps = ceil(T/k);
k = T/numsteps;
[Matrix_solution_Krylov,vectork,evaluateF_1,vec_MVMs_1,vec_inner_product_1,vec_Scalar_multiplication_1,vec_linear_combination_1,vec_err1_L2,vec_err1_max_L2,vec_err1_max,vec_fetch1,vec_store1] = Run_Simulation_Krylov(T,theta,h,numsteps,index_work,k,internalPoints,nu,cases,size_rho,U0,U0_ori,U_ref,CFL,num_k);

n_Jacobianw1 = (vec_MVMs_1)/16;
vec_MVMs_1 = vec_MVMs_1 + 10*evaluateF_1;
Total_Fetch_Krylov = 18*n_Jacobianw1 + 12*evaluateF_1 + 6*vec_linear_combination_1 + 3*vec_Scalar_multiplication_1 + 6*alpha*vec_inner_product_1 + vec_fetch1(:,1) ;
Total_Store_Krylov = 3*n_Jacobianw1 + 12*evaluateF_1 + 3*vec_linear_combination_1 + 3*vec_Scalar_multiplication_1 + vec_store1(:,1);
Total_MOPS_Krylov = 21*n_Jacobianw1 + 12*evaluateF_1 + 9*vec_linear_combination_1 + 6*vec_Scalar_multiplication_1 + 6*alpha*vec_inner_product_1 + vec_fetch1(:,1)  + vec_store1(:,1);

%% =========== 2/ Using Leja methods ===============================================

CFL = CFL_ori; 
index_work = index_work_array(2);
CFL = CFL/(2^(index_work-1));
k = CFL*min([0.5*h,(h^2/(8*nu)),(h/(theta))]);
numsteps = ceil(T/k);
k = T/numsteps;
[Matrix_solution_Leja,evaluateF_2,U_Leja_T,vec_MVMs_2,vec_inner_product_2,vec_Scalar_multiplication_2,vec_linear_combination_2,vec_err2_L2,vec_err2_max_L2,vec_err2_max,vec_fetch2,vec_store2] = Run_Simulation_Leja(T,theta,h,numsteps,index_work,k,internalPoints,nu,cases,size_rho,U0,U0_ori,U_ref,CFL,num_k);

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
[Matrix_solution_Leja_4th2S,evaluateF_3_4th2S,U_Leja_T_4th2S,vec_MVMs_3_4th2S,vec_inner_product_3_4th2S,vec_Scalar_multiplication_3_4th2S,vec_linear_combination_3_4th2S,vec_err3_L2_4th2S,vec_err3_max_L2_4th2S,vec_err3_max_4th2S,vec_fetch3_4th2S,vec_store3_4th2S] = Run_Simulation_Leja_4th_2stages(T,theta,h,numsteps,index_work,k,internalPoints,nu,cases,size_rho,U0,U0_ori,U_ref,CFL,num_k);

n_Jacobianw3_4th2S = (vec_MVMs_3_4th2S)/16;
vec_MVMs_3_4th2S = vec_MVMs_3_4th2S + 10*evaluateF_3_4th2S;
Total_Fetch_Leja_4th2S = 18*n_Jacobianw3_4th2S + 12*evaluateF_3_4th2S + 6*vec_linear_combination_3_4th2S + 3*vec_Scalar_multiplication_3_4th2S + 6*alpha*vec_inner_product_3_4th2S + vec_fetch3_4th2S(:,1);
Total_Store_Leja_4th2S = 3*n_Jacobianw3_4th2S + 12*evaluateF_3_4th2S + 3*vec_linear_combination_3_4th2S + 3*vec_Scalar_multiplication_3_4th2S + + vec_store3_4th2S(:,1) ;
Total_MOPS_Leja_4th2S = 21*n_Jacobianw3_4th2S + 12*evaluateF_3_4th2S + 9*vec_linear_combination_3_4th2S + 6*vec_Scalar_multiplication_3_4th2S + 9*alpha*vec_inner_product_3_4th2S + vec_fetch3_4th2S(:,1) + vec_store3_4th2S(:,1);

%% RK4

num_k_RK4 = 5;
[Matrix_solution_RK4,vectork_rk4,evaluateF_4,vec_MVMs_4,vec_inner_product_4,vec_linear_combination_4,vec_linear_combination5vectors_4,vec_err4_L2,vec_err4_max_L2,vec_err4_max,vec_fetch4,vec_store4] = Run_Simulation_RK4(num_k_RK4,T,h,internalPoints,nu,cases,size_rho,U0,U_ref,upper_bound_RK4,module_eig_max);

vec_MVMs_4 = vec_MVMs_4 + 10*evaluateF_4;
Total_Fetch_RK4 = 12*evaluateF_4 + 6*vec_linear_combination_4 + 15*vec_linear_combination5vectors_4  + vec_fetch4(:,1);
Total_Store_RK4 = 12*evaluateF_4 + 3*vec_linear_combination_4 + 3*vec_linear_combination5vectors_4  + vec_store4(:,1);
Total_MOPS_RK4 = 12*evaluateF_4 + 9*vec_linear_combination_4 + 18*vec_linear_combination5vectors_4  + vec_fetch4(:,1) + vec_store4(:,1);

%% RK2

num_k_RK2 = 5;
[Matrix_solution_RK2,vectork_rk2,evaluateF_5,vec_MVMs_5,vec_inner_product_5,vec_linear_combination_5,vec_linear_combination3vectors_5,vec_err5_L2,vec_err5_max_L2,vec_err5_max,vec_fetch5,vec_store5] = Run_Simulation_RK2(num_k_RK2,T,h,internalPoints,nu,cases,size_rho,U0,U_ref,upper_bound_RK2,module_eig_max);

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

%% Error in discrete L2 norm
h_plot = figure;
loglog(vectork, vec_err1_L2,'-g^','LineWidth',2,'MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',5);
hold on
% loglog(vectork, vec_err2_L2,'-rx','LineWidth',2,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',5);
loglog(vectork,vec_err3_L2,'-mo','LineWidth',2,'MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',5);
loglog(vectork,vec_err3_L2_4th2S,'-.mo','LineWidth',2,'MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',5);
loglog(vectork_rk4, vec_err4_L2,'-s','Color',[0.6350, 0.0780, 0.1840],'LineWidth',2,'MarkerEdgeColor',[0.6350, 0.0780, 0.1840],'MarkerFaceColor',[0.6350, 0.0780, 0.1840],'MarkerSize',5);
loglog(vectork_rk2, vec_err5_L2,'-bd','LineWidth',2,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',5);
loglog(vectork, vectork.^2 - 0.95*vectork.^2,'-k.','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',0.001);
loglog(vectork, vectork.^4 - 0.5*vectork.^4,'--k.','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',0.001);
ylim([-inf 10]);
hold off
set(h_plot,'Units','Inches');
pos = get(h_plot,'Position');
set(h_plot,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h_plot,'filename','-dpdf','-r0')
% legend('Krylov','Chebyshev','Leja','Leja4th2S','RK4','RK2','f(\tau) = c\tau^2','f(\tau) = c\tau^4','location','best');
set(gca,'FontSize',15)
set(gca,'TickLength',[0.03,0])
xlabel('Time step size \tau','FontSize', 17);
ylabel('Errors in discrete L^2 norm','FontSize', 17);
legend('Krylov','Leja','Leja4th2S','RK4','RK2','f(\tau) = c\tau^2','f(\tau) = c\tau^4','location','best','NumColumns',2,'FontSize', 13);

%% Total number of memory operations

h_plot = figure;
loglog( vectork, Total_MOPS_Krylov,'-g^','LineWidth',2,'MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',5);
hold on
% loglog( vectork, Total_MOPS_Chebyshev,'-rx','LineWidth',2,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',5);
loglog( vectork, Total_MOPS_Leja,'-mo','LineWidth',2,'MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',5);
loglog( vectork, Total_MOPS_Leja_4th2S,'-.mo','LineWidth',2,'MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',5);
loglog(vectork_rk4, Total_MOPS_RK4(:,1),'-s','Color',[0.6350, 0.0780, 0.1840],'LineWidth',2,'MarkerEdgeColor',[0.6350, 0.0780, 0.1840],'MarkerFaceColor',[0.6350, 0.0780, 0.1840],'MarkerSize',5);
loglog(vectork_rk2, Total_MOPS_RK2(:,1),'-bd','LineWidth',2,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',5);
hold off
set(h_plot,'Units','Inches');
pos = get(h_plot,'Position');
set(h_plot,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h_plot,'filename','-dpdf','-r0')
set(gca,'FontSize',15)
set(gca,'TickLength',[0.03,0])
xlabel('Time step size \tau','FontSize', 17);
ylabel('Total number of memory operations','FontSize', 17);
legend('Krylov','Leja','Leja4th2S','RK4','RK2','f(\tau) = c\tau^2','f(\tau) = c\tau^4','location','best','FontSize', 13);

%% plot the solution

mesh_x = sparse([ax:h:bx-h]'); %  equidistantly spaced in x-direction
mesh_y = mesh_x;
[X,Y] = meshgrid(mesh_x,mesh_y);
X = X.';
Y = Y.';
%% Original
rho0 = U0(1:size_rho);
u0 = U0(size_rho+1:2*size_rho);
v0 = U0(2*size_rho+1:end);
matrix_rho0 =  reshape(rho0,[dof,dof]);
matrix_u0 =  reshape(u0,[dof,dof]);
matrix_v0 =  reshape(v0,[dof,dof]);

%% numerical solution provided by Leja
rho_Leja = U_Leja_T(1:size_rho);
u_Leja = U_Leja_T(size_rho+1:2*size_rho);
v_Leja = U_Leja_T(2*size_rho+1:end);
matrix_rho_Leja =  reshape(rho_Leja,[dof,dof]);
matrix_u_Leja = reshape(u_Leja,[dof,dof]);
matrix_v_Leja = reshape(v_Leja,[dof,dof]);

%% Reference solution
rho_ref = U_ref(1:size_rho);
u_ref = U_ref(size_rho+1:2*size_rho);
v_ref = U_ref(2*size_rho+1:end);
matrix_rho_ref =  reshape(rho_ref,[dof,dof]);
matrix_u_ref = reshape(u_ref,[dof,dof]);
matrix_v_ref = reshape(v_ref,[dof,dof]);

%% Initial data for rho
h_plot = figure;
s = pcolor(X,Y,matrix_rho0);
xlabel('x');
ylabel('y');
set(gca,'FontSize',12)
set(h_plot,'Units','Inches');
pos = get(h_plot,'Position');
set(h_plot,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h_plot,'filename','-dpdf','-r0')
colorbar
s.MeshStyle = 'none';
s.FaceColor = 'interp';

%% Initial data for u
h_plot = figure;
s = pcolor(X,Y,matrix_u0);
xlabel('x');
ylabel('y');
set(gca,'FontSize',12)
set(h_plot,'Units','Inches');
pos = get(h_plot,'Position');
set(h_plot,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h_plot,'filename','-dpdf','-r0')
% title('T=0');
colorbar
s.MeshStyle = 'none';
s.FaceColor = 'interp';

%% Initial data for v
h_plot = figure;
s = pcolor(X,Y,matrix_v0);
xlabel('x');
ylabel('y');
set(gca,'FontSize',12)
set(h_plot,'Units','Inches');
pos = get(h_plot,'Position');
set(h_plot,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h_plot,'filename','-dpdf','-r0')
colorbar
s.MeshStyle = 'none';
s.FaceColor = 'interp';

%% vorticity for initial data
A_nabla = lap1d_nabla(dof,cases);
Dxx = lap1d(dof,cases);
Iz = speye(dof);
h_plot = figure;
vorticity_u0 = kron(Iz,A_nabla)*v0-kron(A_nabla,Iz)*u0;
matrix_vorticity_u0 =  reshape(vorticity_u0,[dof,dof]);
s = pcolor(X,Y,matrix_vorticity_u0);
xlabel('x');
ylabel('y');
set(gca,'FontSize',12)
set(h_plot,'Units','Inches');
pos = get(h_plot,'Position');
set(h_plot,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h_plot,'filename','-dpdf','-r0')
% title('T=0');
colorbar
s.MeshStyle = 'none';
s.FaceColor = 'interp';

%% reference solution for rho at T 
h_plot = figure;
s = pcolor(X,Y,matrix_rho_ref);
xlabel('x');
ylabel('y');
set(gca,'FontSize',12)
set(h_plot,'Units','Inches');
pos = get(h_plot,'Position');
set(h_plot,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h_plot,'filename','-dpdf','-r0')
colorbar
s.MeshStyle = 'none';
s.FaceColor = 'interp';

%% reference solution for u at T
h_plot = figure;
s = pcolor(X,Y,matrix_u_ref);
xlabel('x');
ylabel('y');
set(gca,'FontSize',12)
set(h_plot,'Units','Inches');
pos = get(h_plot,'Position');
set(h_plot,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h_plot,'filename','-dpdf','-r0')
colorbar
s.MeshStyle = 'none';
s.FaceColor = 'interp';

%% reference solution for v at T
h_plot = figure;
s = pcolor(X,Y,matrix_v_ref);
xlabel('x');
ylabel('y');
set(gca,'FontSize',12)
set(h_plot,'Units','Inches');
pos = get(h_plot,'Position');
set(h_plot,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h_plot,'filename','-dpdf','-r0')
colorbar
s.MeshStyle = 'none';
s.FaceColor = 'interp';

%% vorticity for reference solution
A_nabla = lap1d_nabla(dof,cases);
Dxx = lap1d(dof,cases);
Iz = speye(dof);

h_plot = figure;
vorticity_ref = kron(Iz,A_nabla)*v_ref-kron(A_nabla,Iz)*u_ref;
matrix_vorticity_ref =  reshape(vorticity_ref,[dof,dof]);
s = pcolor(X,Y,matrix_vorticity_ref);
xlabel('x');
ylabel('y');
set(gca,'FontSize',12)
set(h_plot,'Units','Inches');
pos = get(h_plot,'Position');
set(h_plot,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h_plot,'filename','-dpdf','-r0')
colorbar
s.MeshStyle = 'none';
s.FaceColor = 'interp';

%% solution for rho at T provided by Leja
h_plot = figure;
s = pcolor(X,Y,matrix_rho_Leja);
xlabel('x');
ylabel('y');
set(gca,'FontSize',15)
set(h_plot,'Units','Inches');
pos = get(h_plot,'Position');
set(h_plot,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h_plot,'filename','-dpdf','-r0')
colorbar
s.MeshStyle = 'none';
s.FaceColor = 'interp';

%% solution for u at T provided by Leja
h_plot = figure;
s = pcolor(X,Y,matrix_u_Leja);
xlabel('x');
ylabel('y');
set(gca,'FontSize',15)

set(h_plot,'Units','Inches');
pos = get(h_plot,'Position');
set(h_plot,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h_plot,'filename','-dpdf','-r0')
colorbar
s.MeshStyle = 'none';
s.FaceColor = 'interp';

%% solution for v at T provided by Leja
h_plot = figure;
s = pcolor(X,Y,matrix_v_Leja);
xlabel('x');
ylabel('y');
set(gca,'FontSize',15)
set(h_plot,'Units','Inches');
pos = get(h_plot,'Position');
set(h_plot,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h_plot,'filename','-dpdf','-r0')
colorbar
s.MeshStyle = 'none';
s.FaceColor = 'interp';

%% vorticity for solution provided by Leja
h_plot = figure;
vorticity_Leja = kron(Iz,A_nabla)*v_Leja-kron(A_nabla,Iz)*u_Leja;
matrix_vorticity_Leja =  reshape(vorticity_Leja,[dof,dof]);
s = pcolor(X,Y,matrix_vorticity_Leja);
xlabel('x');
ylabel('y');
set(gca,'FontSize',15)
set(h_plot,'Units','Inches');
pos = get(h_plot,'Position');
set(h_plot,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h_plot,'filename','-dpdf','-r0')
colorbar
s.MeshStyle = 'none';
s.FaceColor = 'interp';
%%
% [~,full_eig]  = index_2D(C,CFL,internalPoints,case_initial);
Jacobian = Jacobian_matrix(dof,cases,1e-4,size_rho,U0_ori);
Jacobian1 = Jacobian_matrix(dof,cases,1e-6,size_rho,U0_ori);
full_eig  = eig(full(Jacobian));
full_eig1  = eig(full(Jacobian1));
h_plot = figure;
plot(real(full_eig),imag(full_eig),'*')
hold on
plot(real(full_eig1),imag(full_eig1),'*')
set(h_plot,'Units','Inches');
pos = get(h_plot,'Position');
set(h_plot,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h_plot,'filename','-dpdf','-r0')
xlabel('Re(z)');
ylabel('Im(z)');
set(gca,'FontSize',15)
set(gca,'TickLength',[0.025,0])
xlim([-1.5 0.5])
legend('\nu = 10^{-4}','\nu = 10^{-6}')

%% 
% [~,full_eig]  = index_2D(C,CFL,internalPoints,case_initial);
Jacobian = Jacobian_matrix(dof,cases,1e-4,size_rho,U_ref);
Jacobian1 = Jacobian_matrix(dof,cases,1e-6,size_rho,U_ref);
full_eig  = eig(full(Jacobian));
full_eig1  = eig(full(Jacobian1));
h_plot = figure;
plot(real(full_eig),imag(full_eig),'*')
hold on
plot(real(full_eig1),imag(full_eig1),'*')
set(h_plot,'Units','Inches');
pos = get(h_plot,'Position');
set(h_plot,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h_plot,'filename','-dpdf','-r0')
xlabel('Re(z)');
ylabel('Im(z)');
set(gca,'FontSize',15)
set(gca,'TickLength',[0.025,0])
xlim([-1.5 0.5])
legend('\nu = 10^{-4}','\nu = 10^{-6}')







