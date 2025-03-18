function [z] =  Leja_point(t0,T,n)
    
if (n == 1)
    z = T;
elseif (n==2)
    z(1) = T;
    z(2) = t0;
else 

    z = [];
    z_next = 0;
    z = [z z_next];
    z_next = T;
    z = [z z_next];    
    u = zeros(2,2)
    u(:,1) = 1
    u(:,2) = -[2;-2] %-z(1:2)
    w = conv(u(1,:),u(2,:))
    ind = 1:length(w)
    ind = ind - 1
    ind = flip(ind)
    w_prime = w(1:end-1).*ind(1:end-1)
    critical_points = roots(w_prime)
    fx = f_leja(w,critical_points)
    [~, index_max] = max(fx)
    z_next = critical_points(index_max);
    z = [z z_next];
    
%     z = [];
%     z_next = T;
%     z = [z z_next];
%     z_next = t0;
%     z = [z z_next];    
%     u = zeros(2,2);
%     u(:,1) = 1;
%     u(:,2) = -z(1:2);
%     w = conv(u(1,:),u(2,:));
%     ind = 1:length(w);
%     ind = ind - 1;
%     ind = flip(ind);
%     w_prime = w(1:end-1).*ind(1:end-1);
%     critical_points = roots(w_prime);
%     fx = f_leja(w,critical_points);
%     [~, index_max] = max(fx);
%     z_next = critical_points(index_max);
%     z = [z z_next];
    for i = 4:n
        temp = [1 -z(end)];
        w = conv(w,temp);
        ind = 1:length(w);
        ind = ind - 1;
        ind = flip(ind);
        w_prime = w(1:end-1).*ind(1:end-1);
        critical_points = roots(w_prime);
%         critical_points
        image = imag(critical_points);
        [ind] = find(image);
        critical_points(ind) = [];
        ind_out_domain1 = find(critical_points > T);
        critical_points(ind_out_domain1) = [];
        ind_out_domain2 = find(critical_points < t0);
        critical_points(ind_out_domain2) = [];
%         critical_points
        fx = f_leja(w,critical_points);
        [~, index_max] = max(fx);
        critical_points(index_max);
        if isempty(index_max) == 0 
           z_next = critical_points(index_max);
           z = [z z_next];
        end
        critical_points = [];
        ind_out_domain1 =[];
        ind_out_domain2 =[];
        fx = [];
    end 
    
end
















