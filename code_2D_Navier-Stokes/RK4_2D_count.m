function [ U, MVMs,evaluateF,inner_product,linear_combination,Linear_combination5vectors, fetch, store ] = RK4_2D_count(internalPoints,nu,cases,size_rho,t,U0)
    
    fetch = zeros(2,1);
    store = zeros(2,1);
    MVMs = 0;
    inner_product = 0;
    linear_combination = 0;
    Linear_combination5vectors = 0;
    evaluateF = 0;
    
    dof = internalPoints+1;
        
    for k = 1:length(t)-1
        
        dt = t(k+1)-t(k);
        k1 = vectorF_2D(dof,cases,nu,size_rho,U0);
        k2 = vectorF_2D(dof,cases,nu,size_rho,U0+dt*0.5*k1);
        k3 = vectorF_2D(dof,cases,nu,size_rho,U0+dt*0.5*k2);
        k4 = vectorF_2D(dof,cases,nu,size_rho,U0+dt*k3);
        
        U = U0 + (dt/6)*(k1+2*k2+2*k3+k4);
        U0 = U;
        evaluateF = evaluateF + 4;
        linear_combination = linear_combination  + 3;
        Linear_combination5vectors = Linear_combination5vectors + 1;
        fetch(1) = fetch(1) + 3;
        store(1) = store(1) + 3;
    end
end

