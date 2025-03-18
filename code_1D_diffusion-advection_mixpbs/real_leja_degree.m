function [p,MVMs,inner_product,Scalar_multiplication,linear_combination,fetch,store,Operations] = real_leja_degree(tau, A, u0,a,b, tol, cases)

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
fetch = zeros(2,1);
store = zeros(2,1);
gamma = (b-a)/4.0;
c = 0.5*(a+b);

switch cases
    case 0
        fx = exp((c+gamma*x));
    case 1
        fx = (exp((c+gamma*x))-1)./((c+gamma*x)); %*
        fx(2) = 1;
    case 3
        fx = (exp(1i*tau*(c+gamma*x))-1-1i*tau*(c+gamma*x)-0.5*1i*tau*(c+gamma*x).^2)./((1i*tau*(c+gamma*x)).^3);
        fx(3) = 1;
end

x = transpose(x);
fx = transpose(fx);
[d,Operation0] = divided_dif_sparse(x,fx);
m = 1;
inner_product = 0;
linear_combination = 0;
y = u0;
p = u0*d(m);
MVMs = 0;
Scalar_multiplication = 1;
fetch(1) = fetch(1) + 1;
store(1) = store(1) + 1;

for m = 2:max_leja
    alpha = 1.0/gamma;
    beta  = -c/gamma - x(m-1);
%     p_temp = y; %1 fetching, 1 storing
    y =  y*beta+(alpha*tau)*((A*y)); %4 fetching, 2 storing
    p = p+ d(m)*y; %4 fetching, 2 storing
    MVMs = MVMs + 1;
%     inner_product = inner_product + 1;
    linear_combination = linear_combination + 2;
    if sqrt(sum(abs(y).^2)/length(y))*abs(d(m)) < tol %2 fetching
        break;
    end
%     Operations = Operations + nnz(y) + 2*nnz(A) - length(A) + 2*n + 2*n  + 3*n;
%     fetch(1) = fetch(1) + 1;
%     store(1) = store(1) + 1;

end

end