function [a] = defect_identify_1D(Jacobian,k,eig_max,type_eigenvalue)

eig_max = k*eig_max
[extreigs] = gersh(k*Jacobian);


if (strcmp(type_eigenvalue,'REAL') == 1)
    eigenvalue_min_estimated = floor(extreigs.SR);
    eigenvalue_min_exact = floor(real(eig_max));
    a = abs(eigenvalue_min_exact - eigenvalue_min_estimated);
elseif  (strcmp(type_eigenvalue,'COMPLEX') == 1)
    eigenvalue_min_estimated = floor(extreigs.SI);
    if imag(eig_max) > 0
        eig_max = -eig_max;
    end
    eigenvalue_min_exact = floor(imag(eig_max));
    a = abs(eigenvalue_min_exact - eigenvalue_min_estimated);
    if a==0
        a = 1;
    end
end

end
