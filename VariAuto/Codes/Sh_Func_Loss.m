classdef Sh_Func_Loss < handle
    
    properties (Access = private)
        designVariable
        network

        % Pau addition
        data
        batchSize
        Xb
        Yb
    end
    
    methods (Access = public)
        
        function obj = Sh_Func_Loss(cParams)
            obj.init(cParams)            
        end
        
        function [j,dj] = computeCostAndGradient(obj,order, i)        
            % Versió estocàstica (en la normal treure línia 21)
            obj.updateBatch(order, i);
            j  = obj.computeCost();
            dj = obj.computeGradient();            
        end

        % Ha de desaparèixer
        function [nD, nB, batchSize] = getBatchSize(obj)
            nD = obj.data.Batch_nD;
            nB = obj.data.Batch_nB;
            batchSize = obj.data.batchSize;
        end

        % Ha de desaparèixer
        function [Xtest, Ytest] = getTestData(obj)
            Xtest = obj.data.Xtest;
            Ytest = obj.data.Ytest;
        end
        
    end
    
    methods (Access = private)

        function obj = updateBatch(obj, order, i)
            [X, Y] = obj.data.createMinibatch(order, i);
            obj.Xb = X;
            obj.Yb = Y;
        end

        function [x,y,I] = updateMinibatch(obj,order,i)
            % Està WIP
            
            Xl = obj.data.Xtrain;
            Yl = obj.data.Ytrain;
            I = obj.batchSize;

            cont = 1;
            if i == fix(size(Xl,1)/I)
                plus = mod(size(Xl,1),I);
                x = zeros([I+plus,size(Xl,2)]);
                y = zeros([I+plus,size(Yl,2)]);
            else
                plus = 0;
                x = zeros([I,size(Xl,2)]);
                y = zeros([I,size(Yl,2)]);
            end
            for j = (i-1)*I+1:i*I+plus
                x(cont,:) = Xl(order(j),:);
                y(cont,:) = Yl(order(j),:);
                cont = cont+1;
            end

            obj.Xb = X;
            obj.Yb = Y;

        end
        
        function init(obj,cParams)
            obj.designVariable = cParams.designVariable;
            obj.network        = cParams.network;
            obj.data           = cParams.data;
        end

       function j = computeCost(obj)
           j = obj.network.forwardprop(obj.Xb,obj.Yb);
       end

       function dj = computeGradient(obj)
           dj = obj.network.backprop(obj.Yb);
       end                  
        
    end
    
end