function [ehAv,MVMs,Inner_product,Scalar_multiplication,Linear_combination,fetch,store,Operations] = Krylov_sparse(A,u0,m,n,cases,U0,U_temp,tau,atol,rtol)

fetch = zeros(2,1);
store = zeros(2,1);
Operations = 0;

[Q,H,MVMs,Inner_product,Scalar_multiplication,Linear_combination,fetch_Arnoldi,store_Arnoldi,k] = Arnoldi_sparse(A,u0,m,n,U0,U_temp,tau,atol,rtol,cases);

if (k >= n)
    HH = (H(1:k,1:k));
    V = Q;
    e1 = zeros(k,1);
    e1(1) = 1;
else
    HH = (H(1:k,1:k));
    V = Q;
    e1 = zeros(k,1);
    e1(1) = 1;
end

[QQ,D] = eig(full(HH));

for i=1:size(D,1)
    [D(i,i),Operations2] = varphi(D(i,i),cases);
    Operations = Operations + Operations2;
end


ehAv = real(V*(QQ*(D*(QQ\e1))));
fetch(1) = fetch(1) + fetch_Arnoldi(1) + size(V,2);
store(1) = store(1) + store_Arnoldi(1) + 1;
fetch(2) = fetch(2) + fetch_Arnoldi(2);

