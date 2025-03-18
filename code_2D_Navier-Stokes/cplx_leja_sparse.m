function [Time_Leja1,p,m,MVMs,inner_product,linear_combination,fetch,store,Operations] = cplx_leja_sparse(tau, A,a,b, u0, tol, cases, shift)

if shift == 0
    Time_Leja1 = zeros(3,1);
    tstart = tic;
    % original complex Leja
    max_leja = 110;
    lejaptsarray;
    x = (lejapts(1,:));
    x(1) = -2;
    x(2) = 2;
    x(3) = 0;
    x = x(1:max_leja).';
    sizex = length(x);
    
    Operations = 0;
    gamma = (b-a)/4.0;
    c = (0.5*(a + b));
    fetch = zeros(2,1);
    store = zeros(2,1);
    Operations = Operations + 4;
    
    switch cases %*
        case 0
            fx = exp(1i*tau*(c+gamma*x));
            Operations = Operations + 3*sizex + 3;
        case 1
            %fx = sparse((exp(1i*tau*(c+gamma*x))-1)./(1i*tau*(c+gamma*x)));
            fx = (exp(1i*tau*(c+gamma*x))-1)./(1i*tau*(c+gamma*x));
            fx(3) = 1;
            Operations = Operations + 5*sizex + 6;
        case 3
            fx = (exp(1i*tau*(c+gamma*x))-1-1i*tau*(c+gamma*x)-0.5*1i*tau*(c+gamma*x).^2)./((1i*tau*(c+gamma*x)).^3);
            fx(3) = 1;
    end
    
    fx = transpose(fx);
    Time_Leja1(3) = toc(tstart);
    tstart = tic;
    [d,Operation0] = divided_dif_sparse(x,fx);
%     size(d)
    Time_Leja1(1) = toc(tstart);
    Operations = Operations + Operation0;
    n = size(A,1);
    y = u0;
    p = u0;
    inner_product = 0;
    linear_combination = 0;
    m = 1;
    tstart = tic;
    p = p*d(m);
    Operations = Operations + nnz(u0);
    MVMs = 0;
    fetch(1) = fetch(1) + 6;
    store(1) = store(1) + 6;
    
    for m = 2:max_leja
        alpha = 1.0/gamma;
        beta  = -c/gamma - x(m-1);
        p_temp = y;
        y = (y*beta+alpha*A*p_temp*(-1i));
        p = p + d(m)*y;
        MVMs = MVMs + 9;        
        inner_product = inner_product + 1;
        linear_combination = linear_combination + 2;
        if sqrt(sum(abs(y).^2)/length(y))*abs(d(m)) < tol
            break;
        end
        fetch(1) = fetch(1) + 24;
        store(1) = store(1) + 9;
        fetch(2) = fetch(2) + 9;
        Operations = Operations + 4 + nnz(y) + 2*nnz(A) + 2*n + 2*n + 2 + 3*n;
    end
    
    if m >= max_leja-1
        fprintf('max number of Leja iterations reached\n');
    end
        p = real(p);
    %m = m+1; %due to *
    Time_Leja1(2) = toc(tstart);
    
    
elseif shift == 1
    % by the paper of Peter Kandolf
    
    t = tau;
    v = u0;
    lejaptsarray;
    xi = lejapts(1,:);
    xi = 1i*xi;
    
    N = size(A,1);
    
    A_hermit = t*(A+A')/2;
    nu = 0
    alpha = min(eig(A_hermit))
    
    A_skewhermit = t*(A-A')/2;
    beta = max(abs(eig(A_skewhermit)))
    
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
            %         fx = (exp(1i*tau*(c+gamma*x))-1)./(1i*tau*(c+gamma*x));
            %             fxi(1) = 1;
        case 3
            fx = (exp(1i*tau*(c+gamma*x))-1-1i*tau*(c+gamma*x)-0.5*1i*tau*(c+gamma*x).^2)./((1i*tau*(c+gamma*x)).^3);
            fx(3) = 1;
    end
    xi = transpose(xi);
    fxi = transpose(fxi);
    d = divided_dif(xi,fxi);
    
    p = v;
    m = 1;
    
    p = d(m)*p;
    Identity = eye(N);
    r = (1/gamma)*(t*A-c*Identity)*v;
    max_leja = 90;
    
    for m = 2:2:max_leja
        
            q = (1/gamma)*(t*A-c*Identity)*r;
            
            p = p + real(d(m))*r + d(m+1)*q;
            
            qq = (d(m+1)/real(d(m)))*q;
            addvector = r + qq;
            
            r = (1/gamma)*(t*A-c*Identity)*q + (imag(xi(m)))^2*r;
            
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
    x = sparse(lejapts(1,:));
    x(1) = -2;
    x(2) = 2;
    x(3) = 0;
    %     x(3) = [];
%    x;   
%%
%    N = size(A,1);
%    [spectrum_adv] = max(eig(full(A)));
%    b = spectrum_adv;
%    a = conj(b);
%%   
    gamma = (b-a)/4.0;
    c = 0.5*(a + b);
    
    switch cases
        case 0
            fx = sparse(exp(tau*(c+gamma*x)));
        case 1
            fx = sparse((exp(tau*(c+gamma*x))-1)./(tau*(c+gamma*x)));
        case 3
            fx = ((exp(1i*tau*(c+gamma*x))-1-1i*tau*(c+gamma*x)-0.5*1i*tau*(c+gamma*x).^2)./((1i*tau*(c+gamma*x)).^3));
            fx(3) = 1;
    end
    %     x = transpose(1i*(c+gamma*x));
    x = transpose(x);
    fx = transpose(fx);
    d = divided_dif_sparse(x,fx);
    
    y = u0;
    p = u0;
    m = 1;
    MVMs = 1;
    inner_product = 0;
    
    p = p*d(m);
    max_leja = 90;
    
    for m = 2:max_leja
        alpha = 1.0/gamma;
        beta  = -c/gamma - x(m-1);
        p_temp = y;
        y = y*beta;
        y = y+alpha*A*p_temp;
        p = p + d(m)*y;
        MVMs = MVMs + 3;
        inner_product = inner_product + 1;
        
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









