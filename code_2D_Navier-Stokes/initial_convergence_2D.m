function [rho0,u0,v0] = initial_convergence(ax,bx,h,theta,cases,case_initial)

if cases == 1
    h = (bx-ax)/(N_x+1);
    mesh_node = [ax+h:h:bx-h]'; %  \equidistantly spaced in x-direction
    switch case_initial
        case 1
            u0 = exp(-20*(mesh_node-0.5).^2);
            rho0 = exp(-20*(mesh_node-0.5).^2);
        case 2
            u0 = rand(size(mesh_node));
            rho0 = rand(size(mesh_node));
        case 3
            u0 = normpdf(mesh_node,0.5,1); %normal distribution
            rho0 = normpdf(mesh_node,0.5,1);
        case 4
            u0 = sin(mesh_node*pi); 
            rho0 = sin(mesh_node*pi); 
        case 6 
            u0 = exp(-20*(mesh_node-0.5).^2);
            rho0 = exp(-20*(mesh_node-0.5).^2);
            u0 = u0/norm(u0); %by the paper on Krylov Subspace ... 1997
            rho0 = rho0/norm(rho0);
            size_Che = [14,12,10,6,4];
            
    end
        
%     u0 = u0/norm(u0); %by the paper on Krylov Subspace ... 1997
%     rho0 = rho0/norm(rho0);
    
    U0 = [rho0;u0];
    
    
%     U0 = U0/norm(U0);    
elseif cases == 2
    mesh_x = sparse([ax:h:bx-h]'); %  equidistantly spaced in x-direction
    mesh_y = mesh_x;

    switch case_initial
        case 1
            u0 = exp(-20*(mesh_node-0.5).^2);
            rho0 = exp(-20*(mesh_node-0.5).^2);
        case 2
            u0 = rand(size(mesh_node));
            rho0 = rand(size(mesh_node));
        case 3
            u0 = normpdf(mesh_node,0.5,1); %normal distribution
            rho0 = normpdf(mesh_node,0.5,1);
        case 4
            u0 = sin(mesh_node*pi); 
            rho0 = sin(mesh_node*pi); 
        case 5
            u0 = 0.05*exp(-20*(mesh_node-0.5).^2);
            rho0 = ones(size(mesh_node,1),1);
            size_Che = [6,13,10,6,4];
        case 6 
            u0 = exp(-20*(mesh_node-0.5).^2);
            rho0 = exp(-20*(mesh_node-0.5).^2);
            u0 = u0/norm(u0); %by the paper on Krylov Subspace ... 1997
            rho0 = rho0/norm(rho0);
            size_Che = [14,12,10,6,4];
        case 7
            
            [X,Y] = meshgrid(mesh_x,mesh_y);
            rho0 = arrayfun(@function_U0,X,Y,ones(length(mesh_x)));
            rho0 = rho0(:);
            rho0 = sparse(rho0);
            u0 = arrayfun(@function_U0,X,Y,2*ones(length(mesh_x)));
            u0 = u0(:);
            u0 = sparse(theta*u0);
            v0 = arrayfun(@function_U0,X,Y,3*ones(length(mesh_x)));
            v0 = v0(:);
            v0 = sparse(theta*v0);
        case 8
             
            [X,Y] = meshgrid(mesh_x,mesh_y);
            X = X.';
            Y = Y.';
            
            rho0 = arrayfun(@function_U0,X,Y,4*ones(length(mesh_x)));
            rho0 = rho0(:);
            rho0 = sparse(rho0);
            u0 = arrayfun(@function_U0,X,Y,5*ones(length(mesh_x)));
            u0 = u0(:);
            u0 = sparse(u0);
            v0 = arrayfun(@function_U0,X,Y,6*ones(length(mesh_x)));
            v0 = v0(:);
            v0 = sparse(v0);
         case 9
             
            [X,Y] = meshgrid(mesh_x,mesh_y);
            rho0 = arrayfun(@function_U0,X,Y,7*ones(length(mesh_x)));
            rho0 = rho0(:);
            u0 = arrayfun(@function_U0,X,Y,8*ones(length(mesh_x)));
            u0 = u0(:);
            v0 = arrayfun(@function_U0,X,Y,9*ones(length(mesh_x)));
            v0 = v0(:);
     
    end
    
end
