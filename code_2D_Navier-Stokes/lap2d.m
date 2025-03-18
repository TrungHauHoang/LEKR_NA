function [A,Operations]=lap2d(internalPoints,cases)
% A=lap2d(nx,ny) matrix of ?delta in 2d on a grid
% of nx internal points in x and ny internal points in y
% numbered by row. Uses the function kron of matlab
switch cases
    case 1
        A = 1;
    case 2
        dof = internalPoints + 1;
%         dof
        Dxx = lap1d(dof,cases);
        Iz = speye(size(Dxx));
        A = kron(Iz,Dxx)+kron(Dxx,Iz);
        Operations = 3*size(A,1);
end
end