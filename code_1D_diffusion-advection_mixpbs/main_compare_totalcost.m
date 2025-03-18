
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
CFL = 320;
alpha = 1;
numpoints = 10;

%% Compute initial data

mesh_x = linspace(ax, bx, 159); % Generate 1000 points between 0 and 1
kappa = 33/5120 + (31/5120)*(tanh(20*(mesh_x - 0.8)));

CFL_ori = CFL;
[U0,dof,index_work_array] = initial_1D(ax,bx,internalPoints,cases_boundary,case_initial,CFL);
U0_ori = U0;
sizeU0 = length(U0);
k_ori = 0.1;

%% Create Jacobian matrix

A_nabla = lap1d_nabla(dof,cases_boundary);
minus_A_diff_x = -A_nabla;
A_diffusion = lap1d(dof,cases_boundary);
Jacobian = Jacobian_matrix(minus_A_diff_x,A_diffusion,kappa);
n = size(Jacobian,1);

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
logspace_vector_kry = logspace(log10(0.00464), -13, numpoints);
k = k_ori/2/2/2;
numsteps = ceil(T/k);
scaled_Jacobian = k*Jacobian;
for i = 1:numpoints
    
    vec_MVMs_1 = 0;
    vec_inner_product_1 = 0;
    vec_Scalar_multiplication_1 = 0;
    vec_linear_combination_1 = 0;
    vec_fetch1 = zeros(2,1);
    vec_store1 = zeros(2,1);
    
    U0_krylov = U0_ori;
    Utemp = U0_ori;
    for j =1:numsteps
        U_ref_posterior = expm((j*k)*Jacobian)*U0_ori;
        norm_U0 = norm(U0);
        [U1,MVMs3,inner_product,Scalar_multiplication,linear_combination,fetch3,store3,Operations3] = Krylov_sparse(scaled_Jacobian,U0/norm_U0,140,n,0,U0,Utemp,k,logspace_vector_kry(i),logspace_vector_kry(i));
        U1 = norm_U0*U1;
        Utemp = U0;
        U0 = U1;
        %Calculate the costs
        vec_MVMs_1 = vec_MVMs_1 + MVMs3;
        vec_inner_product_1 = vec_inner_product_1 + inner_product + 1;
        vec_Scalar_multiplication_1 = vec_Scalar_multiplication_1 + Scalar_multiplication + 2;
        vec_linear_combination_1 = vec_linear_combination_1 + linear_combination;
        vec_fetch1(1) = vec_fetch1(1) + fetch3(1) + 1;
        vec_store1(1) = vec_store1(1) + store3(1) + 1;
    end
    Matrix_solution_Krylov(:,i) = U1;
    %% Calculate the error
    error_u = abs(U_ref-U1);
    vec_err1_L2(i,1) = (h*sum(error_u.^2))^0.5
    vec_err1_max(i,1) = norm(U_ref-U1,'inf')
    
    Total_MOPS_Krylov(i,1) = 2*vec_MVMs_1 + 3*vec_linear_combination_1 + 2*vec_Scalar_multiplication_1 + 2*alpha*vec_inner_product_1 + vec_fetch1(1) + vec_store1(1);
    Total_MOPS_Krylov_1(i,1) = 2*vec_MVMs_1 + 3*vec_linear_combination_1 + 2*vec_Scalar_multiplication_1 + 2*alpha*10*vec_inner_product_1 + vec_fetch1(1) + vec_store1(1);
    U0 = U0_ori;
    
end

vec_err1_L2 = vec_err1_L2(vec_err1_L2~=0);
Total_MOPS_Krylov = Total_MOPS_Krylov(Total_MOPS_Krylov~=0);

%% =========== 2/ Using Leja methods ===============================================

k = k_ori/2;
numsteps = ceil(T/k);
tolerance = logspace(-4, -14, 10);

for i=1:numpoints
    vec_MVMs_2 = 0;
    vec_inner_product_2 = 0;
    vec_Scalar_multiplication_2 = 0;
    vec_linear_combination_2 = 0;
    vec_fetch2 = zeros(2,1);
    vec_store2 = zeros(2,1);
    U0_leja = U0_ori;
    if (strcmp(type_eigenvalue,'REAL') == 1)

        spectrum_end = min(real(eig(full(Jacobian))));
        eigenvalue_min = k*spectrum_end;
        eigenvalue_max = 0;

        for j =1:numsteps
            norm_U0 = norm(U0);
            [U1,MVMs2,inner_product,Scalar_multiplication,linear_combination,fetch2,store2,Operations2] = real_leja_degree(k,Jacobian,U0/norm_U0,eigenvalue_min,eigenvalue_max,tolerance(i),0);
            U1 = norm_U0*U1;
            U0 = U1;
            vec_MVMs_2 = vec_MVMs_2 + MVMs2;
            vec_inner_product_2 = vec_inner_product_2 + inner_product + 1;
            vec_Scalar_multiplication_2 = vec_Scalar_multiplication_2 + Scalar_multiplication + 2;
            vec_linear_combination_2 = vec_linear_combination_2 + linear_combination;
            vec_fetch2(1) = vec_fetch2(1) + fetch2(1) + 1;
            vec_store2(1) = vec_store2(1) + store2(1) + 1;
        end
    elseif (strcmp(type_eigenvalue,'COMPLEX') == 1)
        
        spectrum_adv = max(abs(imag(eig(full(Jacobian)))));
        eigenvalue_min = -spectrum_adv;
        eigenvalue_max = spectrum_adv;

        for j =1:numsteps
            norm_U0 = norm(U0);
            [U1,MVMs2,inner_product,Scalar_multiplication,linear_combination,fetch2,store2,Operations2] = cplx_leja_degree(k,Jacobian,U0/norm_U0,eigenvalue_min,eigenvalue_max,tolerance(i),0,0);
            U1 = norm_U0*U1;
            U0 = U1;
            vec_MVMs_2 = vec_MVMs_2 + MVMs2;
            vec_inner_product_2 = vec_inner_product_2 + inner_product + 1;
            vec_Scalar_multiplication_2 = vec_Scalar_multiplication_2 + Scalar_multiplication + 2;
            vec_linear_combination_2 = vec_linear_combination_2 + linear_combination;
            vec_fetch2(1) = vec_fetch2(1) + fetch2(1) + 1;
            vec_store2(1) = vec_store2(1) + store2(1) + 1;
        end
    end
    % Calculate the error
    error_u = abs(U_ref-U1);
    vec_err2_L2(i,1) =  (h*sum(error_u.^2))^0.5
    vec_err2_max(i,1) = norm(U_ref-U1,'inf')
    % Returning to original, to compute the next time step
    U0 = U0_ori;
    Total_MOPS_Leja(i,1) = 2*vec_MVMs_2 + 3*vec_linear_combination_2 + 2*vec_Scalar_multiplication_2 + 2*alpha*vec_inner_product_2 + vec_fetch2(1)  + vec_store2(1);
end

vec_err2_L2 = vec_err2_L2(vec_err2_L2~=0);
Total_MOPS_Leja = Total_MOPS_Leja(Total_MOPS_Leja~=0);
%%

k = k_ori/2/2/2;
numsteps = ceil(T/k);
tolerance = logspace(-4, -13, 10);
for i=1:numpoints
    vec_MVMs_3 = 0;
    vec_inner_product_3 = 0;
    vec_Scalar_multiplication_3 = 0;
    vec_linear_combination_3 = 0;
    vec_fetch3 = zeros(2,1);
    vec_store3 = zeros(2,1);
    U0_leja = U0_ori;
    
    if (strcmp(type_eigenvalue,'REAL') == 1)

        spectrum_end = min(real(eig(full(Jacobian))));
        eigenvalue_min = k*spectrum_end;
        eigenvalue_max = 0;

        for j =1:numsteps
            norm_U0 = norm(U0);
            [U1,MVMs3,inner_product,Scalar_multiplication,linear_combination,fetch3,store3,Operations3] = real_leja_degree(k,Jacobian,U0/norm_U0,eigenvalue_min,eigenvalue_max,tolerance(i),0);
% %            Calculate the costs
            U1 = norm_U0*U1;
            U0 = U1;
            vec_MVMs_3 = vec_MVMs_3+ MVMs3;
            vec_inner_product_3 = vec_inner_product_3 + inner_product + 1;
            vec_Scalar_multiplication_3 = vec_Scalar_multiplication_3 + Scalar_multiplication + 2;
            vec_linear_combination_3 = vec_linear_combination_3 + linear_combination;
            vec_fetch3(1) = vec_fetch3(1) + fetch3(1) + 1;
            vec_store3(1) = vec_store3(1) + store3(1) + 1;
        end
    elseif (strcmp(type_eigenvalue,'COMPLEX') == 1)
        
        spectrum_adv = max(abs(imag(eig(full(Jacobian)))));
        eigenvalue_min = -spectrum_adv;
        eigenvalue_max = spectrum_adv;

        for j =1:numsteps
            norm_U0 = norm(U0);
            [U1,MVMs3,inner_product,Scalar_multiplication,linear_combination,fetch3,store3,Operations3] = cplx_leja_degree(k,Jacobian,U0/norm_U0,eigenvalue_min,eigenvalue_max,tolerance(i),0,0);
            U1 = norm_U0*U1;
            U0 = U1;
            vec_MVMs_3 = vec_MVMs_3+ MVMs3;
            vec_inner_product_3 = vec_inner_product_3 + inner_product + 1;
            vec_Scalar_multiplication_3 = vec_Scalar_multiplication_3 + Scalar_multiplication + 2;
            vec_linear_combination_3 = vec_linear_combination_3 + linear_combination;
            vec_fetch3(1) = vec_fetch3(1) + fetch3(1) + 1;
            vec_store3(1) = vec_store3(1) + store3(1) + 1;
        end
    end
    % Calculate the error
    error_u = abs(U_ref-U1);
    vec_err2_L2_1(i,1) =  (h*sum(error_u.^2))^0.5
    vec_err2_max(i,1) = norm(U_ref-U1,'inf')
    % Returning to original, to compute the next time step
    U0 = U0_ori;
    Total_MOPS_Leja_1(i,1) = 2*vec_MVMs_3 + 3*vec_linear_combination_3 + 2*vec_Scalar_multiplication_3 + 2*alpha*vec_inner_product_3 + vec_fetch3(1) + vec_store3(1);
end

vec_err2_L2_1 = vec_err2_L2_1(vec_err2_L2_1~=0);
Total_MOPS_Leja_1 = Total_MOPS_Leja_1(Total_MOPS_Leja_1~=0);

%% RK4

U0 = U0_ori;
num_k_RK4 = 5;
step_sizerk4 = logspace(log10(upper_bound_RK4/module_eig_max), log10(upper_bound_RK4/(module_eig_max*16)), num_k_RK4);
numsteps_arr4 = ceil(T./step_sizerk4);
step_sizerk4 = T./numsteps_arr4;

for i = 1:num_k_RK4
    vec_fetch4 = zeros(2,1);
    vec_store4 = zeros(2,1);
    
    t = 0:step_sizerk4(i):T;
    [ U1, MVMs1,linear_combination,Linear_combination5vectors, fetch, store, Operations ] = RK4_1D_count(Jacobian,t,U0);

    error_u = abs(U_ref-U1);
    vec_err4_L2(i,1) = (h*sum(error_u.^2))^0.5
    vec_err4_max(i,1) = norm(U_ref-U1,'inf')
    
    vec_fetch4(1) = fetch(1);
    vec_store4(1) = store(1)+3;
    vec_fetch4(2) = fetch(2);
    vec_store4(2) = store(2);
    vec_MVMs_4 = MVMs1;
    vec_linear_combination_4 = linear_combination;
    vec_linear_combination5vectors_4 = Linear_combination5vectors;
    
    Total_MOPS_RK4(i,1) = 2*vec_MVMs_4 + 3*vec_linear_combination_4 + 6*vec_linear_combination5vectors_4  + vec_fetch4(1) + 5*vec_fetch4(2) + vec_store4(1) + 5*vec_store4(2);
end

%% RK2
U0 = U0_ori;
num_k_RK2 = 5;
step_sizerk2 = logspace(log10(upper_bound_RK2/module_eig_max), log10(upper_bound_RK2/(module_eig_max*16)), num_k_RK2);
numsteps_arr2 = ceil(T./step_sizerk2);
step_sizerk2 = T./numsteps_arr2;

for i = 1:num_k_RK2
    
    vec_fetch5 = zeros(2,1);
    vec_store5 = zeros(2,1);
    t = 0:step_sizerk2(i):T;
    [ U1, MVMs1,linear_combination,Linear_combination3vectors, fetch, store, Operations  ] = RK2_1D_count(Jacobian,t,U0);
    error_u = abs(U_ref-U1);
    vec_err5_L2(i,1) = (h*sum(error_u.^2))^0.5
    vec_err5_max(i,1) = norm(U_ref-U1,'inf')
    vec_fetch5(1) = fetch(1);
    vec_store5(1) = store(1)+1;
    vec_fetch5(2) = fetch(2);
    vec_store5(2) = store(2);
    vec_MVMs_5 = MVMs1;
    vec_linear_combination_5 = linear_combination;
    vec_linear_combination3vectors_5 = Linear_combination3vectors;
    Total_MOPS_RK2(i,1) = 2*vec_MVMs_5 + 3*vec_linear_combination_5 + 4*vec_linear_combination3vectors_5  + vec_fetch5(1) + 5*vec_fetch5(2) + vec_store5(1) + 5*vec_store5(2);
    
end

%% =========== compare ====================================================
%% Error in discrete L2 norm
h_plot = figure;
loglog( Total_MOPS_Krylov, vec_err1_L2,'-g^','LineWidth',1.5,'MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',5);
hold on
loglog( Total_MOPS_Krylov_1, vec_err1_L2,'--g^','LineWidth',1.5,'MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',5);
loglog(Total_MOPS_Leja,vec_err2_L2,'--mo','LineWidth',1.5,'MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',5);
loglog(Total_MOPS_Leja_1,vec_err2_L2_1,'-mo','LineWidth',1.5,'MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',5);
% loglog(Total_MOPS_Leja_2,vec_err3_L2_2,'-mo','LineWidth',1.5,'MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',5);
loglog(Total_MOPS_RK4, vec_err4_L2,'-s','Color',[0.6350, 0.0780, 0.1840],'LineWidth',1.5,'MarkerEdgeColor',[0.6350, 0.0780, 0.1840],'MarkerFaceColor',[0.6350, 0.0780, 0.1840],'MarkerSize',5);
loglog(Total_MOPS_RK2, vec_err5_L2,'-bd','LineWidth',1.5,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',5);
ylim([-inf 1e-1]);
hold off
set(h_plot,'Units','Inches');
pos = get(h_plot,'Position');
set(h_plot,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h_plot,'filename','-dpdf','-r0')
set(gca,'FontSize',13)
set(gca,'TickLength',[0.03,0])
xlabel('Total number of memory operations','FontSize', 17);
ylabel('Errors in discrete L^2 norm','FontSize', 17);
xlim([1000 2*1000000])
legend('Krylov \zeta=1','Krylov \zeta=10','Leja \tau=1/20','Leja \tau=1/80','RK4','RK2','location','best','FontSize', 13);

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



