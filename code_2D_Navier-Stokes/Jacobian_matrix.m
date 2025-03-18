
%%========================================================================

function Jacobian = Jacobian_matrix(dof,cases,nu,size_rho,U)

A_nabla = lap1d_nabla(dof,cases);
Dxx = lap1d(dof,cases);
Iz = speye(dof);

% Jacobian = [-kron(Iz,A_nabla)*diag(sparse(U(size_rho+1:2*size_rho)))-kron(A_nabla,Iz)*diag(sparse(U(2*size_rho+1:end))), -kron(Iz,A_nabla)*diag(sparse(U(1:size_rho))), -kron(A_nabla,Iz)*diag(sparse(U(1:size_rho))); ...
%     (1./(sparse(U(1:size_rho).^2))).*diag(kron(Iz,A_nabla)*sparse(U(1:size_rho)))-(1./(sparse(U(1:size_rho)))).*kron(Iz,A_nabla),-(-kron(A_nabla,Iz)*diag(sparse(U(2*size_rho+1:end)))).'+diag(-kron(Iz,A_nabla)*sparse(U(size_rho+1:2*size_rho)))-(-kron(Iz,A_nabla)*diag(sparse(U(size_rho+1:2*size_rho)))).'+nu*(kron(Iz,Dxx)+kron(Dxx,Iz)),diag(-kron(A_nabla,Iz)*sparse(U(size_rho+1:2*size_rho))); ...
%     (1./(sparse(U(1:size_rho).^2))).*diag(kron(A_nabla,Iz)*sparse(U(1:size_rho)))-(1./(sparse(U(1:size_rho)))).*kron(A_nabla,Iz),diag(-kron(Iz,A_nabla)*sparse(U(2*size_rho+1:end))),-(-kron(Iz,A_nabla)*diag(sparse(U(size_rho+1:2*size_rho)))).'+diag(kron(A_nabla,Iz)*sparse(U(2*size_rho+1:end)))-(-kron(A_nabla,Iz)*diag(sparse(U(2*size_rho+1:end)))).'+nu*(kron(Iz,Dxx)+kron(Dxx,Iz))];

Jacobian = [-kron(Iz,A_nabla)*diag(sparse(U(size_rho+1:2*size_rho)))-kron(A_nabla,Iz)*diag(sparse(U(2*size_rho+1:end))), -kron(Iz,A_nabla)*diag(sparse(U(1:size_rho))), -kron(A_nabla,Iz)*diag(sparse(U(1:size_rho))); ...
    (1./(sparse(U(1:size_rho).^2))).*diag(kron(Iz,A_nabla)*sparse(U(1:size_rho)))-(1./(sparse(U(1:size_rho)))).*kron(Iz,A_nabla), -kron(Iz,A_nabla)*diag(sparse(U(size_rho+1:2*size_rho))) - sparse(U(2*size_rho+1:end)).*kron(A_nabla,Iz) + nu*(kron(Iz,Dxx)+kron(Dxx,Iz)),diag(-kron(A_nabla,Iz)*sparse(U(size_rho+1:2*size_rho))); ...
    (1./(sparse(U(1:size_rho).^2))).*diag(kron(A_nabla,Iz)*sparse(U(1:size_rho)))-(1./(sparse(U(1:size_rho)))).*kron(A_nabla,Iz),diag(-kron(Iz,A_nabla)*sparse(U(2*size_rho+1:end))),-kron(A_nabla,Iz)*diag(sparse(U(2*size_rho+1:end))) - sparse(U(size_rho+1:2*size_rho)).*kron(Iz,A_nabla) + nu*(kron(Iz,Dxx)+kron(Dxx,Iz))];

end
