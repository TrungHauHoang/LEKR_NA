function A=lap1d_nabla(n,cases)

if cases == 1 % Homogeneous Dirichlet boundary condition
    h = 1/(n+1); %Because in the case of Homogeneous Dirichlet boundary condtion
    % we don't have one more unknown.
    e = ones(n,1);
    A = spdiags([-e/(2*h) 0*e e/(2*h)],-1:1,n,n);
%     A = spdiags([-e/h e/h],-1:0,n,n);
    
%     A(1,1) = -1/h;
%     A(1,2) = 1/h;
%     A(end,end) = 1/h;
%     A(end,end-1) = -1/h;
    A = full(A);
    
elseif cases == 2 % Periodic boundary condition
    
    h = 1/(n); %Because in the case of periodic boundary condtion
    % we have one more unknown.
    e = ones(n,1);
    A = (spdiags([-e/(2*h) 0*e e/(2*h)],-1:1,n,n));
    A(1,end) = -1/(2*h);
    A(end,1) = 1/(2*h);
    
end

