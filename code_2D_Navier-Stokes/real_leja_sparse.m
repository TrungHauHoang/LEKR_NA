function [Time_Leja1,p,m,MVMs,inner_product,fetch,store,Operations] = real_leja_sparse(tau, A,a,b, u0, tol, cases)
    Time_Leja1 = zeros(3,1);
    tstart = tic;
    max_leja = 110;
    lejaptsarray;
    x = (lejapts(1,:));
    x(1) = -2;
    x(2) = 2;
    x(3) = 0;
    x = x(1:max_leja);
    sizex = length(x);
    n = size(A,1);
    Operations = 0;
    gamma = (b-a)/4.0;
    c = 0.5*(a+b);
    fetch = zeros(2,1);
    store = zeros(2,1);
    Operations = Operations + 4;
    
    switch cases %*
        case 0
            fx = exp((c+gamma*x));
            Operations = Operations + 3*sizex;
        case 1
            fx = (exp((c+gamma*x))-1)./((c+gamma*x)); %*
            fx(2) = 1;
            Operations = Operations + 5*sizex;
        case 3
            fx = (exp(1i*tau*(c+gamma*x))-1-1i*tau*(c+gamma*x)-0.5*1i*tau*(c+gamma*x).^2)./((1i*tau*(c+gamma*x)).^3);
            fx(3) = 1;
    end

    x = transpose(x);
    fx = transpose(fx);
    Time_Leja1(3) = toc(tstart);
    tstart = tic;
    [d,Operation0] = divided_dif_sparse(x,fx);
    Time_Leja1(1) = toc(tstart);
    Operations = Operations + Operation0;
    y = u0;
    
    m = 1;
    tstart = tic;
    inner_product = 0;

    p = u0*d(m);
    Operations = Operations + nnz(u0);
    MVMs = 3;
    
    fetch(1) = fetch(1) + 6;
    store(1) = store(1) + 6;
    
   for m = 2:max_leja
        alpha = 1.0/gamma;
        beta  = -c/gamma - x(m-1);
        p_temp = y; %2 fetching, 2 storing
%         y = y; %2 fetching, 2 storing
        y = y*beta+(alpha*tau)*((A*p_temp)); %4 fetching, 2 storing
        p = p+ d(m)*y; %4 fetching, 2 storing
        MVMs = MVMs + 12;
        inner_product = inner_product + 1;
        if sqrt(sum(abs(y).^2)/length(y))*abs(d(m)) < tol %2 fetching
            break;
        end
%         Operations = Operations + 4 + n + 2*n^2 + 2*n + 2*n + 2 + 3*n;
        Operations = Operations + 4 + nnz(y) + 2*nnz(A) + 2*n + 2*n + 2 + 3*n;
        
        fetch(1) = fetch(1) + 24;
        store(1) = store(1) + 9;
        fetch(2) = fetch(2) + 9;
    end
    Time_Leja1(2) = toc(tstart);
    if m >= max_leja-1
        fprintf('max number of Leja iterations reached\n');
    end
%    m = m + 1;
end

%% Original version
% function [p,m] = real_leja_sparse(tau, A, u0, tol, cases)
% 
%     lejaptsarray;
%     x = sparse(lejapts(1,:));
%     x(1) = -2;
%     x(2) = 2;
%     x(3) = 0;
%     
%     N = size(A,1);
%     
%     
% %     spectrum_end = max(real(eig(A))); 
% %     a = 0;
% %     b = tau*spectrum_end;
% 
%     spectrum_end = min(real(eig(full(A))));
%     a = tau*spectrum_end;
%     b = 0;
%     
%     
%     gamma = (b-a)/4.0;
%     c = 0.5*(a+b);
% %     c
% %     gamma
% 
%     switch cases %*
%         case 0
%             fx = sparse(exp((c+gamma*x)));
%         case 1
%             fx = sparse((exp((c+gamma*x))-1)./((c+gamma*x))); %*
%             fx(2) = 1;
%         case 3
%             fx = sparse((exp(1i*tau*(c+gamma*x))-1-1i*tau*(c+gamma*x)-0.5*1i*tau*(c+gamma*x).^2)./((1i*tau*(c+gamma*x)).^3));
%             fx(3) = 1;
%     end
% 
% %     fx
% %     fx = exp(c+gamma*x);
%     x = transpose(x);
%     fx = transpose(fx);
%     d = divided_dif_sparse(x,fx);
% 
%     y = u0;
%     p = u0;
%     m = 1;
%      %due to *
%     p = p*d(m);
%     max_leja = 90;
% 
%    for m = 2:max_leja
%         alpha = 1.0/gamma;
%         beta  = -c/gamma - x(m-1);
% 
%         p_temp = y;
%         y = y*beta;
%         y = y+alpha*tau*A*p_temp;
%         p = p+ d(m)*y;
% 
%         if sqrt(sum(abs(y).^2)/length(y))*abs(d(m)) < tol
%             break;
%         end
%     end
%     
%     if m >= max_leja-1
%         fprintf('ERROR: max number of Leja iterations reached\n');
%     end
%     m = m + 1;
% end