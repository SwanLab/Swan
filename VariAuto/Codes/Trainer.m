classdef Trainer < handle

    properties (Access = public)
        isDisplayed
        costHist
        optHist
    end
    
    properties (Access = protected) 
       data
       figureCost
       figureOpt
       xIter
       nPlot
    end

    properties (Access = protected)
       optimizationProblem
    end

    methods (Access = public, Static)

        function self = create(varargin)
           switch varargin{2}
               case 'SGD'
                   self = SGD(varargin);
               case 'Fminunc'
                   self = Fminunc(varargin);
               case 'Nesterov'
                   self = Nesterov(varargin);
               case 'RMSProp'
                   self = RMSProp(varargin);
           end
        end
    end

    methods (Access = protected)

        function init(self,s)
            self.optimizationProblem = s{1};
            if length(s) <= 7
                self.isDisplayed  = false;
            else
                self.isDisplayed  = true;
            end
        end

        function [J,g] = costFunction(self,x,Xb,Yb)
            theta   = x;
            self.optimizationProblem.computeCost(theta,Xb,Yb)
            J = self.optimizationProblem.cost;
            g = self.optimizationProblem.gradient;
        end

        function storeValues(self,x,f,state,opt)
            switch state
                case 'init'
                    self.costHist = [0,0,0];
                    self.optHist = [0,0];      
                    self.figureCost = figure;
                    self.figureOpt = figure;
                case 'iter'
                    cV = zeros(1,3);
                    cV(1) = f;
                    cV(2) = self.optimizationProblem.regularization;
                    cV(3) = self.optimizationProblem.loss;
                    self.xIter = [self.xIter, x];
                    self.costHist = [self.costHist;cV];
                    oV = zeros(1,2);
                    oV(1) = opt.gnorm;
                    oV(2) = opt.epsilon;
                    self.optHist = [self.optHist;oV];
            end
        end

        function plotMinimization(self,iter)
            nIter = self.nPlot;
            v = 0:nIter:iter-1;
            if iter > 1
            self.plotCostRegErr(v);
            self.plotEpsOpt(v)
            end
            if self.optimizationProblem.data.nFeatures <= 2 %CUIDADO, Arreglar
                self.optimizationProblem.plotBoundary('contour')
            end
        end  

        function plotCostRegErr (self,v)
            figure(self.figureCost)
            %semilogy(v,self.costHist(2:end,1),'+-r',v,self.costHist(2:end,3),'+-b',v,self.costHist(2:end,2),'+-k')
            plot(v(2:end),self.costHist(2:end,1),'d--b','MarkerFaceColor','b')
            %legend('Fval','Loss','Regularization')
            xlabel('Iterations')
            ylabel('Function Values')
            title('Cost minimization')
            xlim([1,inf])
            drawnow
            hold off
        end

        function plotEpsOpt(self,v)
            figure(self.figureOpt)
%             subplot(2,1,1)
%             plot(v,self.optHist(2:end,1),'+-k')
%             yline(1,'-','Gtol Criteria')
%             xlabel('Iterations')
%             ylabel('Optimalty criteria')
%             title('Gradient norm vs iter')
%             xlim([10,inf])
%             subplot(2,1,2)
            plot(v(2:end),self.optHist(2:end,2),'d--g','MarkerFaceColor','g')
            xlabel('Iterations')
            ylabel('Learning Rate')
            title('Step Size vs Iter')
            xlim([1,inf])
        end
    end
end