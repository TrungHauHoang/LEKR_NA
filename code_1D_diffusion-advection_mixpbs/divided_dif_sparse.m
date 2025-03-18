% function [coef,Operations] =  divided_dif_sparse(x,fx)
%     n = size(x,1);
%     F = sparse(n,n);
%     F(:,1) = fx;
%     Operations = 0;
%     for i=2:n
%         for j=2:i
%             F(i,j) = (F(i,j-1)-F(i-1,j-1))/(x(i)-x(i+1-j));
%             Operations = Operations + 7;
%         end
%     end
%     coef = (diag(F));
% end

function [d,Operations] =  divided_dif_sparse(x,fx)
    d = zeros(length(x),1);
    
    d(1) = fx(1);
    Operations = 0;
    for m=2:length(x)
        d(m) = fx(m);
        for i=1:m-1
            d(m) = (d(m)-d(i))/(x(m)-x(i));
            Operations = Operations + 3;
        end
    end
end