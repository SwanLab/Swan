classdef testingCellParameters < handle
    
    
    methods (Access = public)
        
        function obj = testingCellParameters()
            
            obj.compute()
        end
        
    end
    
    
    methods (Access = private)
        
        function compute
            alpha = sym('alpha');
            rho = sym('rho');
            txi = sym('txi');
            
            mx = sym('mx');
            my = sym('my');
            
            rhoEq = 1 - alpha*mx*my;
            txiEq = 1/2 + (mx - my)/2;
            eq1 = rho - rhoEq;
            eq2 = txi - txiEq;
            
            sol = solve([eq1,eq2],[mx,my]);
            mx = sol.mx;
            my = sol.my;
            
            eqMxMax = mx - 1;
            alphaV = solve(eqMxMax,alpha);
            txiV = solve(alphaV - alpha,txi);
            txiMxMax = matlabFunction(txiV);
            
            eqMyMax = my - 1;
            alphaV = solve(eqMyMax,alpha);
            txiV = solve(alphaV - alpha,txi);
            txiMyMax = matlabFunction(txiV);
            
        end
        
        
    end
    
    
    
end