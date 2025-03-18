function initial = function_U0(x,y,cases) % Defines initial conditions.

switch cases
    case 1
        initial = 0.05*exp(-100*(x-0.5).^2)*exp(-100*(y-0.5).^2)+1;
    case 2
        initial = exp(-20*(x-0.5).^2)*exp(-20*(y-0.5).^2);
        %initial = x+y;
    case 3
        initial = exp(-40*(x-0.5).^2)*exp(-40*(y-0.5).^2);
        %initial = x+y;
    case 4
        initial = 1;
    case 5
        v0 = 0.1;
        Delta = 1/30;
        %         y
        if (0<= y) && (y <= 0.5)
            initial = v0*tanh((y-0.25)/Delta);
        elseif y>0.5
            initial = v0*tanh((0.75-y)/Delta);
        end
    case 6
        delta = 5e-3;
        initial = delta*sin(2*pi*x);
    case 7
        R = 10^-2;
        initial = 1/10;
        if (x^2+y^2) <= R^2
            initial = 1;
        end
    case 8
        initial = 0;
    case 9
        initial = 0;
        
end
