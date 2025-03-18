function [ U, MVMs,Linear_combination,Linear_combination5vectors, fetch, store ] = RK4_2D_count(Jacobian,t,U0)
    
    fetch = zeros(2,1);
    store = zeros(2,1);
    MVMs = 0;
    Linear_combination = 0;
    Linear_combination5vectors = 0;
    
    
    for k = 1:length(t)-1

        dt = t(k+1)-t(k);
        
        [k1, MVMs1] = vectorF_1D(Jacobian,U0);
        [k2, MVMs2] = vectorF_1D(Jacobian,U0+dt*0.5*k1);
        [k3, MVMs3] = vectorF_1D(Jacobian,U0+dt*0.5*k2);
        [k4, MVMs4] = vectorF_1D(Jacobian,U0+dt*k3);
        U = U0 + (dt/6)*(k1+2*k2+2*k3+k4);
        U0 = U;
        MVMs = MVMs + MVMs1 + MVMs2 + MVMs3 + MVMs4;        
        Linear_combination = Linear_combination  + 3;
        Linear_combination5vectors = Linear_combination5vectors + 1;
        fetch(1) = fetch(1) + 1;
        store(1) = store(1) + 1;
    end
    
end

