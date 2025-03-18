function A=lap1d(n,cases)

% lap1d one dimensional finite difference approximation
% A=lap1d(n) computes a (full) matrix finite difference
% approximation of the one dimensional operator (Delta) on the
% domain Omega=(0,1) using n interior points


if cases == 1
    h=1/(n+1);
    e=ones(n,1);
    A=(spdiags([e/h^2 -2/h^2*e e/h^2],-1:1,n,n));
    
elseif cases == 2
    h=1/(n);
    e=ones(n,1);
    A=(spdiags([e/h^2 -2/h^2*e e/h^2],-1:1,n,n));
    A(1,end) = 1/h^2;
    A(end,1) = 1/h^2;
end


