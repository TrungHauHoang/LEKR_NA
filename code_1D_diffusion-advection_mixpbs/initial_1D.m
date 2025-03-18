function [u0,dof,index_work_array] = initial_1D(ax,bx,internalPoints,cases_boundary,case_initial,CFL)

if cases_boundary == 1
    dof = internalPoints;
    h = (bx-ax)/(internalPoints+1);
    mesh_node = [ax+h:h:bx-h]'; %  \equidistantly spaced in x-direction
    switch case_initial
        case 1
            u0 = mesh_node.*(1-mesh_node);
    end
    
    index_work_array = [1 1];
    
end