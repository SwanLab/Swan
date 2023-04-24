classdef Trainer < handle

    properties (Access = public)
        isDisplayed
        costHist
        optHist
        data
    end
    
    properties (Access = protected) 
       figureCost
       figureOpt
       xIter
       nPlot
       CostFunction
    end

    methods (Access = public, Static)

        function obj = create(varargin)
           switch varargin{2}
               case 'SGD'
                   obj = SGD(varargin);
               case 'Fminunc'
                   obj = Fminunc(varargin);
               case 'Nesterov'
                   obj = Nesterov(varargin);
               case 'RMSProp'
                   obj = RMSProp(varargin);
           end
        end
    end

    methods (Access = protected)

        function init(obj,s)
            obj.CostFunction = s{1};
            if length(s) <= 7
                obj.isDisplayed  = false;
            else
                obj.isDisplayed  = true;
            end
        end

        function [J,g] = costFunction(obj,x,Xb,Yb)
            theta   = x;
            obj.CostFunction.computeCost(theta,Xb,Yb)
            J = obj.CostFunction.cost;
            g = obj.CostFunction.gradient;
        end

        function storeValues(obj,x,f,state,opt)
            switch state
                case 'init'
                    obj.costHist = [0,0,0];
                    obj.optHist = [0,0];      
                    obj.figureCost = figure;
                    obj.figureOpt = figure;
                case 'iter'
                    cV = zeros(1,3);
                    cV(1) = f;
                    cV(2) = obj.CostFunction.regularization;
                    cV(3) = obj.CostFunction.loss;
                    obj.xIter = [obj.xIter, x];
                    obj.costHist = [obj.costHist;cV];
                    oV = zeros(1,2);
                    oV(1) = opt.gnorm;
                    oV(2) = opt.epsilon;
                    obj.optHist = [obj.optHist;oV];
            end
        end

%         function plotMinimization(obj,iter)
%             nIter = obj.nPlot;
%             v = 0:nIter:iter-1;
%             if iter > 1
%             obj.plotCostRegErr(v);
%             obj.plotEpsOpt(v)
%             end
%             if obj.optimizationProblem.data.nFeatures <= 2 %CUIDADO, Arreglar
%                 obj.optimizationProblem.plotBoundary('contour')
%             end
%         end  

        function plotCostRegErr (obj,v)
            figure(obj.figureCost)
            %semilogy(v,obj.costHist(2:end,1),'+-r',v,obj.costHist(2:end,3),'+-b',v,obj.costHist(2:end,2),'+-k')
            plot(v(2:end),obj.costHist(2:end,1),'d--b','MarkerFaceColor','b')
            %legend('Fval','Loss','Regularization')
            xlabel('Iterations')
            ylabel('Function Values')
            title('Cost minimization')
            xlim([1,inf])
            drawnow
            hold off
        end

        function plotEpsOpt(obj,v)
            figure(obj.figureOpt)
%             subplot(2,1,1)
%             plot(v,obj.optHist(2:end,1),'+-k')
%             yline(1,'-','Gtol Criteria')
%             xlabel('Iterations')
%             ylabel('Optimalty criteria')
%             title('Gradient norm vs iter')
%             xlim([10,inf])
%             subplot(2,1,2)
            plot(v(2:end),obj.optHist(2:end,2),'d--g','MarkerFaceColor','g')
            xlabel('Iterations')
            ylabel('Learning Rate')
            title('Step Size vs Iter')
            xlim([1,inf])
        end
    end
end