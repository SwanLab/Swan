classdef testingCellParameters < handle
    
    
    properties (Access = private)
        alpha
        rho
        txi
        mx
        my
        txiMxMax
        txiMyMax
        txiMxMin
        txiMyMin        
    end
    
    methods (Access = public)
        
        function obj = testingCellParameters()
            
            obj.compute()
        end
        
    end
    
    
    methods (Access = private)
        
        function compute(obj)
            obj.init()
            obj.computeMxMy();
           % obj.solveMxLowerBound();
           % obj.solveMyLowerBound();           
            obj.solveMxUpperBound();
            obj.solveMyUpperBound();
            obj.plotCellDomain()
        end
        
        function init(obj)
            obj.alpha = sym('alpha');
            obj.rho = sym('rho');
            obj.txi = sym('txi');
        end
        
        function computeMxMy(obj)
            mxS = sym('mx');
            myS = sym('my');            
            rhoEq = 1 - obj.alpha*mxS*myS;
            txiEq = 1/2 + (mxS - myS)/2;
            eq1 = obj.rho - rhoEq;
            eq2 = obj.txi - txiEq;
            sol = solve([eq1,eq2],[mxS,myS]);
            obj.mx = sol.mx;
            obj.my = sol.my;
        end
        
        function solveMxLowerBound(obj)
            eqMxMin = obj.mx - 0;
            obj.txiMxMin = obj.solveEquation(eqMxMin);
        end        
        
        function solveMyLowerBound(obj)
            eqMxMin = obj.my - 0;
            obj.txiMxMin = obj.solveEquation(eqMxMin);
        end
        
        function solveMxUpperBound(obj)
            eqMxMax = obj.mx - 1;
            obj.txiMxMax = obj.solveEquation(eqMxMax);
        end        
        
        function solveMyUpperBound(obj)
            eqMyMax = obj.my - 1;
            obj.txiMyMax = obj.solveEquation(eqMyMax);
        end
        
        function txiFun = solveEquation(obj,eq)
            alphaV = solve(eq,obj.alpha);
            txiV = solve(alphaV - obj.alpha,obj.txi);
            txiFun = matlabFunction(txiV);
        end
        
        function plotCellDomain(obj)
            rhoPlot = linspace(0,1,100);
            alphaPlot = pi/4;
            figure(1)
            plot(rhoPlot,obj.txiMxMax(alphaPlot,rhoPlot),'r')
            hold on
            plot(rhoPlot,obj.txiMyMax(alphaPlot,rhoPlot),'r')            
        end
        
    end
    
    
    
end