classdef Trainer < handle

    properties (Access = public)
        isDisplayed
        designVariable
        objectiveFunction
    end
    
    properties (Access = protected) 
        xIter
        nPlot
        %Xtrain
        %Ytrain
        %Xtest
        %Ytest
    end

    properties (Access = private)
        figureOpt
        figureCost
        costHist
        optHist
    end


    methods (Access = public, Static)

        function obj = create(cParams)
           switch cParams.type
               case 'SGD'
                   obj = SGD(cParams);
               case 'Fminunc'
                   obj = Fminunc(cParams);
               case 'Nesterov'
                   obj = Nesterov(cParams);
               case 'RMSProp'
                   obj = RMSProp(cParams);
           end
        end
    end

    methods (Access = protected)

        function init(obj,cParams)
            obj.objectiveFunction = cParams.costFunc;
            obj.designVariable = cParams.designVariable;
            obj.isDisplayed  = false;
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
                    cV(2) = obj.objectiveFunction.regularization;
                    cV(3) = obj.objectiveFunction.loss;
                    obj.xIter = [obj.xIter, x];
                    obj.costHist = [obj.costHist;cV];
                    oV = zeros(1,2);
                    oV(1) = opt.gnorm;
                    oV(2) = opt.epsilon;
                    obj.optHist = [obj.optHist;oV];
            end
        end

        function plotCostRegErr(obj,v)
            figure(obj.figureCost)
            plot(v(2:end),obj.costHist(2:end,1),'d--b','MarkerFaceColor','b')
            xlabel('Iterations')
            ylabel('Function Values')
            title('Cost minimization')
            xlim([1,inf])
            drawnow
            hold off
        end

        function plotEpsOpt(obj,v)
            figure(obj.figureOpt)
            plot(v(2:end),obj.optHist(2:end,2),'d--g','MarkerFaceColor','g')
            xlabel('Iterations')
            ylabel('Learning Rate')
            title('Step Size vs Iter')
            xlim([1,inf])
        end
    end
end