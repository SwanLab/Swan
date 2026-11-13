function [alpha,cw] = createHybridFunction(type,initDerivZero,x0,coeffs)
    syms x
    switch type
        case 'Quad'
            if initDerivZero
                a = -1/(x0^2 - 2*x0); b = 0;
                m = 2*a*x0;
            else
                a = coeffs(1);
                b = a*x0^2 - 2*a*x0 + 1;
                m = 2*a*x0 + b;
            end
            alpha = piecewise((x<x0),(a*x^2 + b*x),(x>=x0),(m*x-m+1));
            cw = double(4*int(sqrt(alpha),0,1));
        case 'Cubic'
            if initDerivZero
                a = coeffs(1); c = 0;
                b = (2*a*x0^3 - 3*a*x0^2 + 1)/(-x0^2 + 2*x0); 
                m = 3*a*x0^2 + 2*b*x0;
            else
                a = coeffs(1); b = coeffs(2);
                c = 2*a*x0^3 + (b-3*a)*x0^2 - 2*b*x0 + 1;
                m = 3*a*x0^2 + 2*b*x0 + c;
            end
            alpha = piecewise((x<x0),(a*x^3 + b*x^2 + c*x),(x>=x0),(m*x-m+1));
            cw = double(4*int(sqrt(alpha),0,1));
    end

end
