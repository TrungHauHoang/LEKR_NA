function [Q,H,MVMs,inner_product,Scalar_multiplication,linear_combination,fetch,store,j] = Arnoldi_sparse(dof,cases,nu,size_rho,k,U,q1,m,cases_varphi,atol)

n = size_rho*3;
MVMs = 0;
Q = zeros(n,m);
Q(:,1) = q1;
fetch = 3;
store = 3;
H = zeros(min(m+1,m),n);
inner_product = 0;
linear_combination = 0;
Scalar_multiplication = 0;
for j=1:m
    z = k*Jacobian_matrix_actvec(dof,cases,nu,size_rho,U,Q(:,j));
    MVMs = MVMs + 16;
    for i=1:j
        H(i,j) = Q(:,i)'*z;
        inner_product = inner_product + 1;
        z = (z - H(i,j)*Q(:,i));
        linear_combination = linear_combination + 1;
    end
    if j < n
        H(j+1,j) = norm(z);
        inner_product = inner_product + 1;
        if H(j+1,j) == 0, return, end
        Q(:,j+1) = z/H(j+1,j);
        Scalar_multiplication = Scalar_multiplication + 1;
    end
    
    HH = H(1:j,1:j);
    [QQ,D] = eig(full(HH));
    for i=1:size(D,1)
        D(i,i) = varphi(D(i,i),cases_varphi+1);
    end
    temp = QQ*D*QQ^-1;
    if  real(H(j+1,j)*temp(j,1)) < atol
        break;
    end
    
end
Q = Q(:,1:j);

