function [Q,H,MVMs,Inner_product,Scalar_multiplication,Linear_combination,k] = Arnoldi_sparse(A,q1,m,n,atol,cases)

Q = zeros(n,m);
Q(:,1) = q1;
H = zeros(min(m+1,m),n);

MVMs = 0;
Inner_product = 0;
Scalar_multiplication = 0;
Linear_combination = 0;

for k=1:m
    z = A*Q(:,k);
    MVMs = MVMs + 1;
    for i=1:k
        H(i,k) = Q(:,i)'*z;
        Inner_product = Inner_product + 1;
        z = (z - H(i,k)*Q(:,i));
        Linear_combination = Linear_combination + 1;
    end
    if k < n
        H(k+1,k) = norm(z);
        Inner_product = Inner_product + 1;
        if H(k+1,k) == 0, return, end
        Q(:,k+1) = z/H(k+1,k);
        Scalar_multiplication = Scalar_multiplication + 1;
    end
        
    %% Stopping criteria Jitse Niesen and Will M. Wright.
        HH = H(1:k,1:k);
        [QQ,D] = eig(full(HH));
        for i=1:size(D,1)
            D(i,i) = varphi(D(i,i),cases+1);
        end
        temp = QQ*D*QQ^-1;
        if  real(H(k+1,k)*temp(k,1)) < atol
            break;
        end
        
        
end
Q = Q(:,1:k);


