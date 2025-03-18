function value = f_leja(w,x)
    n = length(w);
    ind = 1:length(w);
    ind = ind - 1;
    ind = flip(ind);
    value = 0;
    for i=1:n
        temp = w(i)*x.^(ind(i));
        value = value + temp;
    end 
    value = abs(value);
end