classdef ImplicitEquationSolver < handle
    
    properties (Access = public)
        xsol
    end
    
    properties (Access = private)
        x0
        functionToSolve
    end
    
    methods (Access = public)
        
        function obj =  ImplicitEquationSolver(cParams)
            obj.init(cParams)
        end
        
        function r = solve(obj)
            F = obj.functionToSolve;
            [rub,rlb] = obj.findRbounds(F);
            r = fzero(F,[rlb,rub]);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.x0 = cParams.x0;
            obj.functionToSolve = cParams.functionToSolve;
        end        
        
        function [rub,rlb] = findRbounds(obj,F)
            r0 = 1;
            F0 = F(r0);
            eps = 10^(-6);
            if F0 >= 0
                r1 = 1 + eps;
                F1 = F(r1);
                while F1 >= 0
                    rnew = obj.newPointBySecant(r0,r1,F0,F1);
                    r0 = r1;
                    F0 = F1;
                    r1 = rnew;
                    F1 = F(r1);
                end
                rub = r1;
                rlb = r0;
            else
                r0 = 1;
                r1 = 1 - eps;
                F1 = F(r1);
                while F1 <= 0
                    rnew = obj.newPointBySecant(r0,r1,F0,F1);
                    r0 = r1;
                    F0 = F1;
                    r1 = max(1e-12,rnew);
                    F1 = F(r1);
                end
                rub = r0;
                rlb = r1;
            end
        end
        
    end
    
    methods (Access = private, Static)
        
        function x2 = newPointBySecant(x0,x1,f0,f1)
            x2 = x1 - (x1-x0)/(f1 - f0)*f1;
        end
        
    end
    
end