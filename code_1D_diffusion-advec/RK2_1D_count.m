function [ U, MVMs,linear_combination,Linear_combination3vectors, fetch, store ] = RK2_2D_count(Jacobian,t,U0)
    
    fetch = zeros(2,1);
    store = zeros(2,1);
    MVMs = 0; 
    linear_combination = 0;
    Linear_combination3vectors = 0;
        
    for k = 1:length(t)-1
        dt = t(k+1)-t(k);
        [k1, MVMs1] = vectorF_1D(Jacobian,U0);
        [k2, MVMs2] = vectorF_1D(Jacobian,U0+dt*k1);
        U = U0 + 0.5*dt*(k1+k2);
        U0 = U;
        MVMs = MVMs + MVMs1 + MVMs2;  
        linear_combination = linear_combination  + 1;
        Linear_combination3vectors = Linear_combination3vectors + 1;
        fetch(1) = fetch(1)  + 1;
        store(1) = store(1)  + 1;
    end
end