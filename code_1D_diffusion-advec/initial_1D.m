function [u0,dof,index_work_array,num_k,tolerence_kry,tolerence_leja] = initial_1D(ax,bx,internalPoints,cases_boundary,case_initial,CFL)

if cases_boundary == 1
    dof = internalPoints;
    h = (bx-ax)/(internalPoints+1);
    mesh_node = [ax+h:h:bx-h]'; %  \equidistantly spaced in x-direction
    switch case_initial
        case 1
            u0 = mesh_node.*(1-mesh_node);
    end
    
    if CFL == 80 || CFL == 160 || CFL == 640 || CFL == 2.5 || CFL == 320 || CFL == 20 %it run for comparison of total computational cost
        index_work_array = [1 1 1];
        num_k = 0;
        tolerence_kry = 0;
        tolerence_leja = 0;
    end
    
    if CFL == 1280
        index_work_array = [1 3];
        num_k = 14;
        tolerence_kry = [2*1e-2 2*1e-2 1e-2 1e-2 1e-2/5 1e-3 1e-3/4 1e-18 1e-5 1e-5 1e-5 1e-5/5 1e-6 1e-6/5];
        tolerence_leja = [0  0 1e-5 1e-5/2 1e-5/2 1e-4/5 1e-5/2 0  0 1e-8 1e-8 1e-8 1e-7/8 1e-7/10]; 
        
    end
    
    if CFL == 40
        index_work_array = [1 1];
        num_k = 10;
        tolerence_kry = [1e-2/8 1e-3/2 1e-3/3.5 1e-3/6 1e-4/2 1e-5/5 1e-5/7 1e-6/2 1e-6/4 1e-6/6];
        tolerence_leja = [1e-4 1e-3/3 1e-4/2 1e-4/3 1e-4/3 1e-6 2*1e-7 1e-7/1.05 1e-7/5 1e-7/4]; 
        
    end
end