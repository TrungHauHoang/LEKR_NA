function [ehAv,MVMs,inner_product,Scalar_multiplication,linear_combination,fetch,store] = Krylov_sparse(dof,cases_bc,nu,size_rho,k,U,F,m,cases_varphi,atol)

fetch = 0;
store = 0;
n = size(U,1);
[Q,H,MVMs,inner_product,Scalar_multiplication,linear_combination,fetch_Arnoldi,store_Arnoldi,size_Hermit] = Arnoldi_sparse(dof,cases_bc,nu,size_rho,k,U,F,m,cases_varphi,atol);
if (m >= n)
    HH = (H(1:m,1:m));
    V = Q;
    e1 = zeros(m,1);
    e1(1) = 1;
else
    HH = H(1:size_Hermit,1:size_Hermit);
    V = Q;
    e1 = zeros(size_Hermit,1);
    e1(1) = 1;
end
[QQ,D] = eig(full(HH));
for i=1:size(D,1)
    D(i,i) = varphi(D(i,i),cases_varphi);
end
ehAv = real(V*(QQ*(D*(QQ\e1))));
fetch = fetch + fetch_Arnoldi(1) + size(V,2);
store = store + store_Arnoldi(1) + 1;
