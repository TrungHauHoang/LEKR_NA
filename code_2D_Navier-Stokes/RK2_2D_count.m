function [ U, MVMs,evaluateF,inner_product,linear_combination,Linear_combination3vectors, fetch, store ] = RK2_2D_count(internalPoints,nu,cases,size_rho,t,U0)
    
    fetch = zeros(2,1);
    store = zeros(2,1);
    MVMs = 0; 
    inner_product = 0;
    evaluateF = 0;
    linear_combination = 0;
    Linear_combination3vectors = 0;
    size_U = length(U0);
    dof = internalPoints+1;
    
    
    for k = 1:length(t)-1
        dt = t(k+1)-t(k);
        
        k1 = vectorF_2D(dof,cases,nu,size_rho,U0);
        k2 = vectorF_2D(dof,cases,nu,size_rho,U0+dt*k1);
        
        U = U0 + 0.5*dt*(k1+k2);
        U0 = U;
        evaluateF = evaluateF + 2;
        linear_combination = linear_combination + 1; 
        Linear_combination3vectors = Linear_combination3vectors + 1;
        fetch(1) = fetch(1) + 3;
        store(1) = store(1) + 3;
                       
    end
    
end