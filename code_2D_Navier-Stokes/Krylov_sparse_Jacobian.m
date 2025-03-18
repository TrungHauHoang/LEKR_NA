function [Time_Krylov,ehAv,MVMs,inner_product,Scalar_multiplication,linear_combination,fetch,store,Operations] = Krylov_sparse_Jacobian(k,A,u0,m,cases)

fetch = zeros(2,1);
store = zeros(2,1);
Operations = 0;
n = size(A,1);

Time_Krylov = zeros(3,1);
tstart = tic;
[Q,H,MVMs,inner_product,Scalar_multiplication,linear_combination,fetch_Arnoldi,store_Arnoldi,Operations1] = Arnoldi_sparse(k*A,u0,m);
Time_Krylov(1) = toc(tstart);
tstart = tic;

if (m >= n)
    HH = (H(1:m,1:m));
    V = Q;
    e1 = zeros(m,1);
    e1(1) = 1;
else
    HH = (H(1:m,1:m));
    V = Q;
    V(:,end) = [];
    e1 = zeros(m,1);
    e1(1) = 1;
end

% HH

[QQ,D] = eig(full(HH));
Time_Krylov(2) = toc(tstart);
tstart = tic;

for i=1:size(D,1)
    [D(i,i),Operations2] = varphi(D(i,i),cases);
    Operations = Operations + Operations2;
end

ehAv = real(V*(QQ*(D*(QQ\e1))));
Time_Krylov(3) = toc(tstart);
fetch(1) = fetch(1) + fetch_Arnoldi(1) + size(V,2);
store(1) = store(1) + store_Arnoldi(1) + 1;
fetch(2) = 0;
store(2) = 0;
%%
% line 5: n*n 
% line 31: QQ\e1: (2/3)*m^3 + 2*m^2
% line 31: D*(QQ\e1): 2*m^2
% line 31: QQ*(D*(QQ\e1)): 2*m^2
% line 31: V*(QQ*(D*(QQ\e1))): size(V,1)*2*m
Operations = Operations + Operations1 + nnz(A) + (2/3)*m^3 + 2*m^2 + 2*m^2 + 2*m^2 + size(V,1)*2*m;