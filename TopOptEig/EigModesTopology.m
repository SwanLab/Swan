classdef EigModesTopology < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        lambda
    end
    
    methods (Access = public)
        
        function obj = EigModesTopology(cParams)
            obj.init(cParams)
            % createLHS?
        end
        
        function fx = provideFunction(obj)
            obj.computeEigenModesAndValues();
            obj.lambda = obj.computeLambda();
            gamma = obj.designVariable.getFirstEigenMode();
            fx = gamma-obj.lambda(1);
        end
        
        function grad = provideDerivative(obj,eigNum)
            obj.reorderModes(obj.lambda,obj.V,obj.D);
            Belem =  obj.bendingMatComputer.elementalBendingMatrix;
            x = obj.designVariable.getColumnArea();
            nElem = obj.mesh.nelem;
            eigV1 = obj.D(1,1);
            eigV2 = obj.D(2,2);
            difEigs = abs(eigV2-eigV1);
            if difEigs > 1
                dfdx = obj.computeSimpleEig(Belem,x);
            else
                dfdx = obj.computeDoubleEig(Belem,x);
            end
            dfdx(1,nElem+1) = 1;
            dfdx(2,nElem+1) = 1;
            grad = dfdx(eigNum,:);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            
        end
        
        function computeEigenModesAndValues(obj)
            obj.stiffnessMatComputer.compute();
            obj.bendingMatComputer.compute();
            [Kfree,free]  = obj.stiffnessMatComputer.provideFreeStiffnessMatrix();
            obj.freeNodes = free;
            Bfree  = obj.bendingMatComputer.provideFreeBendingMatrix();
            obj.computeEigenFunctionAndValues(Bfree,Kfree);
        end
        
        function l = computeLambda(obj)
            l = sort(diag(obj.D));
        end
    end
    
end