classdef SimpleShapeOptimizationSolver < handle

    properties (Access = public)
        tV
        JV
        incXvalues
        designVariable
    end
    
    properties (Access = private)
        topOptProblem
        cost
        plotter
        TOL
        maxIter
    end
    
    methods (Access = public)
        
        function obj = SimpleShapeOptimizationSolver(cParams)
            obj.init(cParams);
        end
        
        function solve(obj)
            xNew = obj.computeInitialValue();
            xOld = xNew;
            x = xOld;
            incX = obj.computeIncX(xOld,xNew);
            iter = 1;
            while ~obj.hasConverged(iter,incX)
                [J,dJ] = obj.cost.computeValueAndGradient(xNew);
                t = obj.computeLineSearch(x,dJ);
                x = obj.computeGradientStep(x,dJ,t);
                x = obj.computeProjection(x);
                [xOld,xNew] = obj.updateXnewXold(xNew,x);
                incX = obj.computeIncX(xOld,xNew);
                obj.plotCostAndLineSearch(iter,J,t,incX);
                iter = iter + 1;
            end
        end
        
    end
    
    methods (Access = private, Static)
       
        function x = computeGradientStep(x,dJ,t)
            x = x - t*dJ;
        end
        
        function x = computeProjection(x)
           x = max(min(x,1),0);
        end
        
        function [xOld,xNew] = updateXnewXold(xNew,xNewNew)
            xOld = xNew;
            xNew = xNewNew;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.TOL            = cParams.TOL;
            obj.maxIter        = cParams.maxIter;
            obj.topOptProblem  = cParams.topOptProblem;
            obj.createDesignVariable();
            obj.createCost();
            obj.createPlotter();
        end
        
        function createDesignVariable(obj)
            obj.designVariable = obj.topOptProblem.designVariable;
        end
        
        function createCost(obj)
            s.topOpt = obj.topOptProblem;
            s.designVariable = obj.designVariable;
            obj.cost = CostComplianceVolume(s);
        end
        
        function createPlotter(obj)
            s.designVariable = obj.designVariable;
            obj.plotter = PlotterDensity(s);
        end
        
        function x0 = computeInitialValue(obj)
            x0 = obj.designVariable.value;
        end
        
        function incX = computeIncX(obj,xOld,xNew)
           incX  = obj.computeNorm(xOld - xNew);
           xNorm = obj.computeNorm(xOld);
           incX = incX/xNorm;
        end
       
        function itHas = hasConverged(obj,iter,incX)
            if iter == 1
                itHas = false;
            else
                itHas = iter >= obj.maxIter || incX < obj.TOL;
            end
        end
        
        function t = computeLineSearch(obj,x,dJ)
            tC = 100;
            incT = 1;
            tA = obj.computeAdimensionalLineSearch(x,dJ);
            t = max(tA,incT*tC);
        end
        
        function t0 = computeAdimensionalLineSearch(obj,x,dJ)
            nX  = obj.computeNorm(x);
            ndJ = obj.computeNorm(dJ);
            t0 = ndJ/nX;
        end
        
        function plotCostAndLineSearch(obj,iter,J,t,incX)
            obj.JV(iter) = J;
            obj.tV(iter) = t;
            obj.incXvalues(iter) = incX;
            obj.plotter.plot(obj.JV,obj.tV,obj.incXvalues);
        end
        
        function n = computeNorm(obj,x)
            sc = obj.designVariable.scalarProduct;
            s  = sc.computeSP_M(x,x);
            n  = sqrt(s);
        end
        
    end
    
end

