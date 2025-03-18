function [ehAv,MVMs,Inner_product,Scalar_multiplication,Linear_combination,fetch,store] = Krylov_sparse(A,u0,m,n,cases,atol)

fetch = 0;
store = 0;

[Q,H,MVMs,Inner_product,Scalar_multiplication,Linear_combination,k] = Arnoldi_sparse(A,u0,m,n,atol,cases);

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
    D(i,i) = varphi(D(i,i),cases);
end

ehAv = real(V*(QQ*(D*(QQ\e1))));
fetch(1) = fetch(1) + size(V,2);
store(1) = store(1) + 1;
