function [p,MVMs,inner_product,Scalar_multiplication,linear_combination,fetch,store,Operations] = cplx_leja_degree(tau, A, u0,a,b, tol, cases, shift)

if shift == 0
    % original complex Leja
    max_leja = 200;
    list_leja_points;
    x = -x;
    x(1) = 0;
    x(2) = 2;
    x(3) = -2;
    x = x(1:max_leja);
    sizex = length(x);
    n = size(A,1);
    Operations = 0;
    gamma = (b-a)/4.0;
    c = 0.5*(a + b);
    fetch = zeros(2,1);
    store = zeros(2,1);
    
    switch cases %*
        case 0
            fx = exp(1i*tau*(c+gamma*x));
        case 1
            fx = (exp(1i*tau*(c+gamma*x))-1)./(1i*tau*(c+gamma*x));
            fx(3) = 1;
        case 3
            fx = (exp(1i*tau*(c+gamma*x))-1-1i*tau*(c+gamma*x)-0.5*1i*tau*(c+gamma*x).^2)./((1i*tau*(c+gamma*x)).^3);
            fx(3) = 1;
    end
    
    x = transpose(x);
    fx = transpose(fx);
    d = divided_dif_sparse(x,fx);
    y = u0;
    m = 1;
    inner_product = 0;
    linear_combination = 0;
    p = u0*d(m);
    MVMs = 0;
    Scalar_multiplication = 1;
    fetch(1) = fetch(1) + 1;
    store(1) = store(1) + 1;
    
    for m = 2:max_leja
        alpha = 1.0/gamma;
        beta  = -c/gamma - x(m-1);
        y = y*beta+alpha*A*y*(-1i);
        p = p + d(m)*y;
        MVMs = MVMs + 1;
        linear_combination = linear_combination + 2;
        if sqrt(sum(abs(y).^2)/length(y))*abs(d(m)) < tol %2 fetching
            break;
        end
        
    end
    
    
elseif shift == 1
    % by the paper of Peter Kandolf
    max_leja = 200;
    list_leja_points;
    x = -x;
    x(1) = 0;
    x(2) = 2;
    x(3) = -2;
    x = x(1:max_leja);
    
    fetch = zeros(2,1);
    store = zeros(2,1);
    Operations = 0;
    t = tau;
    v = u0;
%     lejaptsarray;
%     xi = lejapts(1,:);
    xi = 1i*x;
    
    N = size(A,1);
    
    A_hermit = t*(A+A')/2;
    nu = 0;
    alpha = min(eig(A_hermit));
    
    A_skewhermit = t*(A-A')/2;
    beta = max(abs(eig(A_skewhermit)));
    
    a = (nu-alpha)/2;
    gamma = 0.5*sqrt((a^(2/3) + beta^(2/3))*abs(a^(4/3) - beta^(4/3)));
    c = (alpha + nu)/2;
    %     c = -0.5;
    %     c = 0;
    switch cases
        case 0
            fxi = (exp(c + gamma*xi));
        case 1
            fxi = (exp(c + gamma*xi) - 1)./(c + gamma*xi);
        case 3
            fx = (exp(1i*tau*(c+gamma*x))-1-1i*tau*(c+gamma*x)-0.5*1i*tau*(c+gamma*x).^2)./((1i*tau*(c+gamma*x)).^3);
            fx(3) = 1;
    end
    xi = transpose(xi);
    fxi = transpose(fxi);
    d = divided_dif_sparse(xi,fxi);
    inner_product = 0;
    linear_combination = 0;
    MVMs = 0;
    Scalar_multiplication = 1;
    fetch(1) = fetch(1) + 1;
    store(1) = store(1) + 1;
    
    p = v;
    m = 1;
    
    p = d(m)*p;
    Identity = eye(N);
    r = (1/gamma)*(t*A-c*Identity)*v;
    
    for m = 2:2:max_leja
        
            q = (1/gamma)*(t*A-c*Identity)*r;
            
            p = p + real(d(m))*r + d(m+1)*q;
            
            qq = (d(m+1)/real(d(m)))*q;
            addvector = r + qq;
            
            r = (1/gamma)*(t*A-c*Identity)*q + (imag(xi(m)))^2*r;
            MVMs = MVMs + 1;
            linear_combination = linear_combination + 2;
            
            if sqrt(sum(abs(addvector).^2)/length(addvector))*abs(real(d(m))) < tol
                break;
            end
        
        if m >= max_leja-1
            fprintf('ERROR: max number of Leja iterations reached\n');
        end
        %     p = real(p);
    end
    m = m+1; %due to *
    
elseif shift == 2
    % my idea
    lejaptsarray;
    x = lejapts(1,:);
    x(1) = -2;
    x(2) = 2;
    x(3) = 0;
    %     x(3) = [];
    x;
    
    N = size(A,1);
    
    spectrum_adv = max(abs(imag(eig(A))));
    a = -spectrum_adv;
    b = spectrum_adv;
    
    gamma = (b-a)/4.0;
    c = 0.5*(a + b);
    
    A_hermit = (A+A')/2;
    nu = 0;
    alpha = min(eig(A_hermit));
    
    shift = -32;    
    
    switch cases
        case 0
            fx = exp(tau*c + 1i*tau*gamma*x);
        case 1
            fx = (exp(tau*shift + 1i*tau*gamma*x)-exp(tau*shift))./(1i*tau*gamma*x);
            
            fx(3) = 1;
        case 3
            fx = (exp(1i*tau*(c+gamma*x))-1-1i*tau*(c+gamma*x)-0.5*1i*tau*(c+gamma*x).^2)./((1i*tau*(c+gamma*x)).^3);
            fx(3) = 1;
    end
    %     x = transpose(1i*(c+gamma*x));
    x = transpose(x);
    fx = transpose(fx);
    d = divided_dif(x,fx);
    
    y = u0;
    p = u0;
    m = 1;
    
    p = p*d(m);
    max_leja = 90;
    
    for m = 2:max_leja
        alpha = 1.0/gamma;
        beta  = -c/gamma - x(m-1);
        p_temp = y;
        y = y*beta;
        y = y+alpha*A*p_temp*(-1i);
        p = p + d(m)*y;
        if sqrt(sum(abs(y).^2)/length(y))*abs(d(m)) < tol
            break;
        end
        
    end
    
    
    if m >= max_leja-1
        fprintf('ERROR: max number of Leja iterations reached\n');
    end
    p = exp(-tau*shift)*p;
    
    m = m+1; %due to *
elseif shift == 3
    % by the paper Exponential integrator
    lejaptsarray;
    x = lejapts(1,:);
    x(1) = -2;
    x(2) = 2;
    x(3) = 0;
    %     x(3) = [];
    x;
    
    N = size(A,1);
    
    [spectrum_adv, index] = max(eig(full(A)));
    b = spectrum_adv;
    a = conj(b);
    
    gamma = (b-a)/4.0;
    c = 0.5*(a + b);
    
    switch cases
        case 0
            fx = exp(tau*(c+gamma*x));
        case 1
            fx = (exp(tau*(c+gamma*x))-1)./(tau*(c+gamma*x));
        case 3
            fx = (exp(1i*tau*(c+gamma*x))-1-1i*tau*(c+gamma*x)-0.5*1i*tau*(c+gamma*x).^2)./((1i*tau*(c+gamma*x)).^3);
            fx(3) = 1;
    end
    %     x = transpose(1i*(c+gamma*x));
    x = transpose(x);
    fx = transpose(fx);
    d = divided_dif(x,fx);
    
    y = u0;
    p = u0;
    m = 1;
    
    p = p*d(m);
    max_leja = 90;
    
    for m = 2:max_leja
        alpha = 1.0/gamma;
        beta  = -c/gamma - x(m-1);
        p_temp = y;
        y = y*beta;
        y = y+alpha*A*p_temp;
        p = p + d(m)*y;
        if sqrt(sum(abs(y).^2)/length(y))*abs(d(m)) < tol
            break;
        end
        
    end
    
    if m >= max_leja-1
        fprintf('ERROR: max number of Leja iterations reached\n');
    end
    %     p = real(p);

    m = m+1; %due to *
end









