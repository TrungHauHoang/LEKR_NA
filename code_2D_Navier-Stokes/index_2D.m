function [U_ref,full_eig] = index_2D(C,CFL,internalPoints,case_initial,theta,h)

switch case_initial
    case 7
        %temporary, need to change
        switch internalPoints
            case 19
                if C == 20 && CFL == 320
                    a = load('saveU_ref40_real');
                    U_ref = a.U_ref_real_20;
                    full_eig = a.full_eig_real_20;
                elseif C == 1.8 && CFL == 80
                    a = load('saveU_ref40_complex');
                    U_ref = a.U_ref_complex_20;
                    full_eig = a.full_eig_complex_20;
                end
            case 39
                if C == 20 && CFL == 320
                    a = load('saveU_ref40_real');
                    U_ref = a.U_ref;
                    full_eig = a.full_eig;
                elseif C == 1.875 && CFL == 80
                    a = load('saveU_ref40_complex');
                    U_ref = a.U_ref;
                    full_eig = a.full_eig;
                end
            case 79
                if C == 20 && CFL == 640
                    a = load('saveU_ref40_real');
                    U_ref = a.U_ref_real_80;
                    full_eig = a.full_eig_real_80;
                elseif C == 1.9 && CFL == 80
                    a = load('saveU_ref40_complex');
                    U_ref = a.U_ref_complex_80;
                    full_eig = a.full_eig_complex_80;
                end
            case 159
                if C == 20 && CFL == 320
                    a = load('saveU_ref40_real');
                    U_ref = a.U_ref_real_160;
                    %             full_eig = a.full_eig_real_80; %out of memory
                elseif C == 1.9 && CFL == 80
                    a = load('saveU_ref40_complex');
                    U_ref = a.U_ref_complex_160;
                    %             full_eig = a.full_eig_complex_80; %out of memory
                end
        end
    case 8
        switch internalPoints
            case 39
                if C == 20 && CFL == 320
                    a = load('saveU_ref40_real');
                    U_ref = a.U_ref_real_40_ShearFlow;
                    full_eig = a.full_eig_real_40_ShearFlow;
                elseif C == 1.8 && CFL == 80
                    a = load('saveU_ref40_complex');
                    U_ref = a.U_ref_complex_40_ShearFlow;
                    full_eig = a.full_eig_complex_40_ShearFlow;
                end
            case 79
                if C == 20 && CFL == 640
                    a = load('saveU_ref40_real');
                    U_ref = a.U_ref_real_80_ShearFlow;
                    full_eig = a.full_eig_real_80_ShearFlow;
                elseif C == 1.8 && CFL == 80
                    a = load('saveU_ref40_complex');
                    U_ref = a.U_ref_complex_80_ShearFlow;
                    full_eig = a.full_eig_complex_80_ShearFlow;
                end
            case 159
                if C == 20 && CFL == 1280
                    a = load('saveU_ref40_real');
                    U_ref = a.U_ref_real_160_ShearFlow;
                    %             full_eig = a.full_eig_real_80; %out of memory
                elseif C == 1.8 && CFL == 160
                    a = load('saveU_ref40_complex');
                    U_ref = a.U_ref_complex_160_ShearFlow;
                    %             full_eig = a.full_eig_complex_80; %out of memory
                elseif C == 12 && CFL == 640
                    a = load('saveU_ref40_real');
                    U_ref = a.U_ref_real_160_ShearFlow0;
                elseif (C == 0.096) && (CFL == 20) || (C == 0.096) && (CFL == 40) || (C == 0.096) && (CFL == 80) || (C == 0.096) && (CFL == 320)
                    a = load('saveU_ref40_complex');
                    U_ref = a.U_ref_complex_160_ShearFlow1;
                elseif C == 0.192 && CFL == 320
                    a = load('saveU_ref40_complex');
                    U_ref = a.U_ref_complex_160_ShearFlow2;
                elseif  C == 1e-6*((1/theta)*(1/h)) && CFL == 320
                    a = load('saveU_ref40_complex');
                    U_ref = a.U_ref_complex_160_ShearFlow3; %reference solution computed by RK4
                    %                     U_ref = a.U_ref_complex_160_ShearFlow4;   %reference solution computed by Krylov
                    %                     U_ref = a.U_ref_complex_160_ShearFlow5; %reference solution computed by Leja
                elseif  C == 1e-3*((1/theta)*(1/h)) && CFL == 20
                    a = load('saveU_ref40_complex');
                    U_ref = a.U_ref_complex_160_ShearFlow6; 
                end
                
        end
        
end
U_ref = full(U_ref);

end

