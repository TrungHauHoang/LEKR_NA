function [p,MVMs,inner_product,Scalar_multiplication,linear_combination,fetch,store] = cplx_leja_sparse_degree(dof,case_bcs,nu,size_rho,tau, U,a,b, u0, m, cases, shift,tol,w3,w2,w1)

if shift == 0
    max_leja = m;
    list_leja_points;
    x = x(1:max_leja).';
    gamma = (b-a)/4.0;
    c = (0.5*(a + b));
    fetch = zeros(2,1);
    store = zeros(2,1);
    
    switch cases
        case 0
            fx = exp(1i*tau*(c+gamma*x));
        case 1
            fx = (exp(1i*tau*(c+gamma*x))-1)./(1i*tau*(c+gamma*x));
            fx(3) = 1;
        case 3
            fx = (exp(1i*tau*(c+gamma*x))-1-1i*tau*(c+gamma*x)-0.5*(1i*tau*(c+gamma*x)).^2)./((1i*tau*(c+gamma*x)).^3);
            fx(3) = 1/6;
    end
    
    fx = transpose(fx);
    [d] = divided_dif_sparse(x,fx);
    linear_combination = 0;
    m = 1;
    p = u0*d(m);
    MVMs = 0;
    Scalar_multiplication = 1;
    inner_product = 0;
    for m = 2:max_leja
        alpha = 1.0/gamma;
        beta  = -c/gamma - x(m-1);
        u0 = (u0*beta+alpha*Jacobian_matrix_actvec(dof,case_bcs,nu,size_rho,U,u0)*(-1i));
        p = p + d(m)*u0;
        MVMs = MVMs + 16;
        linear_combination = linear_combination + 2;
        
        if sqrt(sum(abs(u0).^2)/length(u0))*abs(d(m)) < tol 
            break;
        end
        
    end
%     m
    p = real(p);
    
elseif shift == 1
    max_leja = m;
    list_leja_points;
    x = x(1:max_leja).';
    gamma = (b-a)/4.0;
    c = (0.5*(a + b));
    fetch = zeros(2,1);
    store = zeros(2,1);
    
    switch cases
        case 0
            fx = exp(1i*tau*(c+gamma*x));
        case 1
            fx = (exp(1i*tau*(c+gamma*x))-1)./(1i*tau*(c+gamma*x));
            fx(3) = 1;
        case 3
            fx = (exp(1i*tau*(c+gamma*x))-1-1i*tau*(c+gamma*x)-0.5*(1i*tau*(c+gamma*x)).^2)./((1i*tau*(c+gamma*x)).^3);
            fx(3) = 1/6;
    end
    
    fx = transpose(fx);
    [d] = divided_dif_sparse(x,fx);
    linear_combination = 0;
    m = 1;
    p = u0*d(m);
    MVMs = 0;
    inner_product = 0;
    Scalar_multiplication = 1;
    
    for m = 2:max_leja
        alpha = 1.0/gamma;
        beta  = -c/gamma - x(m-1);
        u0 = (u0*beta+alpha*Jacobian_matrix_actvec_add(dof,case_bcs,nu,size_rho,U,u0,w3,w2,w1)*(-1i));
        p = p + d(m)*u0;
        MVMs = MVMs + 16;
        linear_combination = linear_combination + 2;
        fetch(1) = fetch(1) + 6;
        if sqrt(sum(abs(u0).^2)/length(u0))*abs(d(m)) < tol 
            break;
        end
        
    end
%     m
    p = real(p);
    
end









