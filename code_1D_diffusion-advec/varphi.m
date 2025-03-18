function [f,Operations]=varphi(lambda, cases)

switch cases
    case 0
        f = exp(lambda);
        Operations = 1;
    case 1
        if abs(lambda) < 1e-7
            f = 1+lambda/2+lambda/6;
            Operations = 5;
        else
            f = (exp(lambda)-1)/lambda;
            Operations = 4;
        end
    case 3
        if abs(lambda) < 1e-7
            f = 1/6+lambda/24+lambda/120;
            Operations = 6;
        else
            f =  (exp(lambda)-1-lambda-(lambda^2)/2)/(lambda^3);
            Operations = 11;
        end
end

end