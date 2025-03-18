function A=lap1d_nabla(n,cases)

% A=lap1d(n) computes a (full) matrix mid point
% approximation of the one dimensional operator (\nabla) on the
% domain Omega=(0,1) using n interior points, 

if cases == 1 % Homogeneous Dirichlet boundary condition
    h = 1/(n+1); %Because in the case of Homogeneous Dirichlet boundary condtion
    % we don't have one more unknown.
    e = ones(n,1);
    A = (spdiags([-e/(2*h) 0*e e/(2*h)],-1:1,n,n));
elseif cases == 2 % Periodic boundary condition
    
    h = 1/(n); %Because in the case of periodic boundary condtion
    % we have one more unknown.
    e = ones(n,1);
    A = (spdiags([-e/(2*h) 0*e e/(2*h)],-1:1,n,n));
    A(1,end) = -1/(2*h);
    A(end,1) = 1/(2*h);
    
end

