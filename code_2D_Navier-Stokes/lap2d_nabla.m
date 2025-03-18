function A=lap2d_nabla(internalPoints,cases)

% A=lap1d(n) computes a (full) matrix mid point
% approximation of the one dimensional operator (\nabla) on the
% domain Omega=(0,1) using n interior points, 

if cases == 1 % Homogeneous Dirichlet boundary condition
    h = 1/(n+1); %Because in the case of Homogeneous Dirichlet boundary condtion
    % we don't have one more unknown.
    e = ones(n,1);
    A = (spdiags([-e/(2*h) 0*e e/(2*h)],-1:1,n,n));
elseif cases == 2 % Periodic boundary condition

    dof = internalPoints+1; 
    A_nabla = lap1d_nabla(dof,cases);
    Iz = speye(dof);
    A = (kron(Iz,A_nabla) + kron(A_nabla,Iz));
    
end

