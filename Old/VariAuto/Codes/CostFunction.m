classdef CostFunction < handle
        
    properties (Access = private)
        network
        regularization
        loss
        lambda
        designVariable
    end

    methods (Access = public)

        function obj = CostFunction(cParams)
            obj.init(cParams);
            obj.createLossFunctional(cParams);
            obj.createRegularizationFunctional();
        end

        function [nD, nB, batchSize] = getBatchSize(obj)
            [nD, nB, batchSize] = obj.loss.getBatchSize();
        end

        function [Xtest, Ytest] = getTestData(obj)
            [Xtest, Ytest] = obj.loss.getTestData();
        end
        
        function [j,dj] = computeCost(obj,theta,order,i)
            % Must create mini batch Xb and Yb out of (order, i)
            %[Xb, Yb] = network.createMiniBatch();
            
            obj.designVariable.thetavec = theta;
            [c,dc] = obj.loss.computeCostAndGradient(order, i); 
            [r,dr] = obj.regularization.computeCostAndGradient();
            l = obj.lambda;
            j = c + l*r;  
            dj = dc + l*dr;
        end 

        function h = getOutput(obj,X)
            h = obj.network.computeLastH(X);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
           obj.network = cParams.network;
           obj.lambda  = cParams.lambda;
           obj.designVariable = cParams.designVariable;
        end

        function createLossFunctional(obj, cParams)
            s.network        = obj.network;
            s.designVariable = obj.designVariable;
            s.data           = cParams.data;
            l = Sh_Func_Loss(s);
            obj.loss = l;
        end

        function createRegularizationFunctional(obj)
            s.designVariable = obj.designVariable;
            r = Sh_Func_L2norm(s);
            obj.regularization = r;
        end
    end   

end

