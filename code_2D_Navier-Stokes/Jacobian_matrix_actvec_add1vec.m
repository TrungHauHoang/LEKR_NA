

function [x] = Jacobian_matrix_actvec_add1vec(dof,cases,nu,size_rho,U,F,w1)

A_nabla = lap1d_nabla(dof,cases);
Dxx = lap1d(dof,cases);
Iz = speye(dof);

x = [full(-kron(Iz,A_nabla)*((sparse(U(size_rho+1:2*size_rho))).*F(1:size_rho))-kron(A_nabla,Iz)*((sparse(U(2*size_rho+1:end))).*F(1:size_rho)) - kron(Iz,A_nabla)*((sparse(U(1:size_rho))).*F(size_rho+1:2*size_rho)) - kron(A_nabla,Iz)*((sparse(U(1:size_rho))).*F(2*size_rho+1:3*size_rho))+ F(end)*w1(1:size_rho));...
    full((1./sparse(U(1:size_rho).^2)).*(kron(Iz,A_nabla)*sparse(U(1:size_rho))).*F(1:size_rho) -(1./sparse(U(1:size_rho))).*(kron(Iz,A_nabla)*F(1:size_rho)) + ((kron(A_nabla,Iz)*diag(sparse(U(2*size_rho+1:end)))).')*F(size_rho+1:2*size_rho) - (kron(Iz,A_nabla)*sparse(U(size_rho+1:2*size_rho))).*F(size_rho+1:2*size_rho)+((kron(Iz,A_nabla)*diag(sparse(U(size_rho+1:2*size_rho)))).')*F(size_rho+1:2*size_rho)+nu*(kron(Iz,Dxx)+kron(Dxx,Iz))*F(size_rho+1:2*size_rho) - (kron(A_nabla,Iz)*sparse(U(size_rho+1:2*size_rho))).*F(2*size_rho+1:3*size_rho) +F(end)*w1(size_rho+1:2*size_rho));...
    full((1./sparse(U(1:size_rho).^2)).*(kron(A_nabla,Iz)*sparse(U(1:size_rho))).*F(1:size_rho) -(1./sparse(U(1:size_rho))).*(kron(A_nabla,Iz)*F(1:size_rho)) - (kron(Iz,A_nabla)*sparse(U(2*size_rho+1:end))).*F(size_rho+1:2*size_rho) +(kron(Iz,A_nabla)*diag(sparse(U(size_rho+1:2*size_rho)))).'*F(2*size_rho+1:3*size_rho)-(kron(A_nabla,Iz)*sparse(U(2*size_rho+1:end))).*F(2*size_rho+1:3*size_rho)+(kron(A_nabla,Iz)*diag(sparse(U(2*size_rho+1:end)))).'*F(2*size_rho+1:3*size_rho)+nu*(kron(Iz,Dxx)+kron(Dxx,Iz))*F(2*size_rho+1:3*size_rho) + F(end)*w1(2*size_rho+1:3*size_rho));...
    full(0)];



