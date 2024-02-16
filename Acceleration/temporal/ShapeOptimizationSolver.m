classdef ShapeOptimizationSolver < handle

    properties (Access = public)
        tV
        JV
        betaV
        incJV
        norm_dJ
        designVariable
    end

    properties (Access = private)
        topOpt
        cost
        plotter
        TOL
        maxIter
        momentumParameter
        momentumParams
        isAccepted
        tau
    end

    methods (Access = public)

        function obj = ShapeOptimizationSolver(cParams)
            obj.init(cParams);
        end

        function solveAdaptative(obj)
            xNew = obj.computeInitialValue();
            xOld = xNew;
            incX = obj.computeIncX(xOld,xNew);
            iter = 1;
            [J,dJ] = obj.cost.computeValueAndGradient(xNew);
            Jold = J;
            incJ = abs(J - Jold);
            t = obj.computeLineSearch(xNew,dJ);
            while ~obj.hasConverged(iter,incJ)
                obj.isAccepted = false;
                % [J,dJ] = obj.cost.computeValueAndGradient(xNew);
                beta = obj.computeBeta(iter);
                y = obj.addMomentumTerm(xOld,xNew,beta);
                y = obj.computeProjection(y);
                while ~obj.isAccepted         
                    [~,dJ] = obj.cost.computeValueAndGradient(y);
                    x = obj.computeGradientStep(y,dJ,t);
                    x = obj.computeProjection(x);
                    [J,~] = obj.cost.computeValueAndGradient(x);
                    % incX = obj.computeIncX(xOld,x);
                    obj.isAccepted = J <= Jold;% && incX < 0.1;
                    if obj.isAccepted
                        t = 1.5*t;
                    elseif ~obj.isAccepted
                        t = t/2;
                        % x = y;
                    end
                    if t < 1e-12
                        error('Failed to converge')
                    end
                end
                [xOld,xNew] = obj.updateXnewXold(xNew,x);
                % incX = obj.computeIncX(xOld,xNew);
                incJ = abs(J - Jold);
                obj.plotCostAndLineSearch(iter,J,t,beta,incJ);
                Jold = J;
                iter = iter + 1;
            end
        end

        function solve(obj)
            xNew = obj.computeInitialValue();
            xOld = xNew;
            [J,dJ] = obj.cost.computeValueAndGradient(xNew);
            J_old = J;
            incJ = 0;
            iter = 1;
            normdJ = norm(dJ);
            while ~obj.hasConverged(iter,incJ)             
                beta = obj.computeBeta(iter);
                x = obj.addMomentumTerm(xOld,xNew,beta);
                t = obj.computeLineSearch(x,dJ);
                x = obj.computeGradientStep(x,dJ,t);
                x = obj.computeProjection(x);
                [xOld,xNew] = obj.updateXnewXold(xNew,x);
                obj.plotCostAndLineSearch(iter,J,t,beta,incJ,normdJ);
                [J,dJ] = obj.cost.computeValueAndGradient(xNew); % This might be needed to move...
                incJ  = abs(J - J_old);
                normdJ = norm(dJ);
                J_old = J;
                iter = iter + 1;
            end
        end

        function solveNesterov(obj)
            xNew = obj.computeInitialValue();
            xOld = xNew;
            [J,dJ] = obj.cost.computeValueAndGradient(xNew);
            J_old = J;
            incJ = 0;
            iter = 1;
            normdJ = norm(dJ);
            while ~obj.hasConverged(iter,incJ)             
                beta = obj.computeBeta(iter);
                y = obj.addMomentumTerm(xOld,xNew,beta);
                [~,dJ] = obj.cost.computeValueAndGradient(y);
                t = obj.computeLineSearch(y,dJ);
                x = obj.computeGradientStep(y,dJ,t);
                x = obj.computeProjection(x);
                [xOld,xNew] = obj.updateXnewXold(xNew,x);
                [J,dJ] = obj.cost.computeValueAndGradient(x);
                obj.plotCostAndLineSearch(iter,J,t,beta,incJ,normdJ);
                incJ = abs(J - J_old);
                normdJ = norm(dJ);
                J_old = J;
                iter = iter + 1;
            end
        end

    end

    methods (Access = private, Static)

        function xNewNew = addMomentumTerm(xOld,xNew,beta)
            xNewNew = xNew + beta*(xNew - xOld);
        end

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

        function beta = computeBeta(obj,iter)
            s.iter = iter;
            beta = obj.momentumParameter.computeValue(s);
        end

        function init(obj,cParams)
            obj.momentumParams = cParams.momentumParams;
            obj.TOL = cParams.TOL;
            obj.maxIter = cParams.maxIter;
            obj.tau = cParams.tau;
            obj.createSettings();
            obj.createDesignVariable();
            obj.createCost();
            % obj.createPlotter();
            obj.createMomentumParameter();
        end

        function createMomentumParameter(obj)
            s = obj.momentumParams;
            obj.momentumParameter = MomentumParameter.create(s);
        end

        function createSettings(obj)
            settings = Settings('Example1');
            translator = SettingsTranslator();
            translator.translate(settings);
            fileName = translator.fileName;
            settingsTopOpt = SettingsTopOptProblem(fileName);
            obj.topOpt = TopOpt_Problem(settingsTopOpt);
        end

        function createDesignVariable(obj)
            obj.designVariable = obj.topOpt.designVariable;
        end

        function createCost(obj)
            s.topOpt = obj.topOpt;
            s.designVariable = obj.designVariable;
            obj.cost = CostComplianceVolume(s);
        end

        function createPlotter(obj)
            s.designVariable = obj.designVariable;
            obj.plotter = PlotterDensity(s);
        end

        function x0 = computeInitialValue(obj)
            x0 = obj.designVariable.fun.fValues;
        end

        function incX = computeIncX(obj,xOld,xNew)
            incX  = obj.computeNorm(xOld - xNew);
            xNorm = obj.computeNorm(xOld);
            incX = incX/xNorm;
        end

        function itHas = hasConverged(obj,iter,norm_dJ)
            if iter == 1
                itHas = false;
            else
                itHas = iter >= obj.maxIter || norm_dJ < obj.TOL;
            end
        end

        function t = computeLineSearch(obj,x,dJ)
            tC = obj.tau;
            incT = 1;
            tA = obj.computeAdimensionalLineSearch(x,dJ);
            t = max(tA,incT*tC);
        end

        function t0 = computeAdimensionalLineSearch(obj,x,dJ)
            nX  = obj.computeNorm(x);
            ndJ = obj.computeNorm(dJ);
            t0 = ndJ/nX;
        end

        function plotCostAndLineSearch(obj,iter,J,t,beta,incJ,norm_dJ)
            obj.JV(iter) = J;
            obj.tV(iter) = t;
            obj.betaV(iter) = beta;
            obj.incJV(iter) = incJ;
            obj.norm_dJ(iter) = norm_dJ;
            % obj.plotter.plot(obj.JV,obj.tV,obj.betaV,obj.incXvalues);
        end

        function n = computeNorm(obj,x)
            % sc = obj.designVariable.scalarProduct;
            % s  = sc.computeSP_M(x,x);
            s  = dot(x,x);
            n  = sqrt(s);
        end

    end

end

