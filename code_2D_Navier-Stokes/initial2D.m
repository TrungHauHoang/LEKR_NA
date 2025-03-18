function [U0,size_rho,dof,index_work_array] = initial2D(ax,bx,internalPoints,h,theta,cases,case_initial,C,CFL)

dof = internalPoints + 1;
mesh_x = ([ax:h:bx-h]'); %  equidistantly spaced in x-direction
mesh_y = mesh_x;

[X,Y] = meshgrid(mesh_x,mesh_y);
X = X.';
Y = Y.';
rho0 = arrayfun(@function_U0,X,Y,4*ones(length(mesh_x)));
rho0 = rho0(:);
u0 = arrayfun(@function_U0,X,Y,5*ones(length(mesh_x)));
u0 = u0(:);
v0 = arrayfun(@function_U0,X,Y,6*ones(length(mesh_x)));
v0 = v0(:);
size_rho = length(rho0);
U0 = [rho0;u0;v0];

if theta == 0.1 ||  theta == 1/6
        index_work_array = [4 1 1];
elseif  theta == 1/12
    if  C == 4.8000e-04  && CFL == 320
        index_work_array = [4 1 1];
    elseif  C == 1e-6*((1/theta)*(1/h)) && CFL == 320
        index_work_array = [4 1 1];
    end
end




