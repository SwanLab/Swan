classdef FE_Interpolation < handle
    
    properties (Access = public)
        shape
        deriv        
        finiteElement
    end
    
    properties (Access = public)
        order
        
        ndime
        nnode
        
        pos_nodes
        isoDv
        
        iteration
        cases
        selectcases
        main_loop
        extra_cases
    end
    
    properties (Access = private)
        type
    end
    
    methods (Access = public)
        
        function obj = FE_Interpolation(cParams)
            obj.init(cParams);
            obj.computeParams();
            obj.computeCases();
        end
        
        function computeShapeDeriv(obj,xGauss)
            obj.computeShapes(xGauss);
            obj.computeShapeDerivatives(xGauss);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.type  = cParams.mesh.type;
            obj.order = cParams.order;
            obj.finiteElement = cParams.finiteElement;
        end
        
        function computeParams(obj)
            obj.ndime = 2;
            obj.nnode = 3;
            obj.pos_nodes = [0 0; 1 0; 0 1];
            obj.isoDv = 0.5;
            obj.main_loop = [3 3];
            obj.extra_cases = [];
        end
        
        function computeShapes(obj,posgp)
            ngaus = size(posgp,2);
%             ndim = 2;
            ndim = 1;
            x = posgp(1,:,:);
            y = posgp(2,:,:);
            fEl = obj.finiteElement;
            obj.shape = zeros(length(fEl.shapeFunctions),ngaus,ndim);

            for s = 1:length(fEl.shapeFunctions)
                obj.shape(s,:,:) = fEl.shapeFunctions{s}(x',y');
            end
        end
        
        function computeShapeDerivatives(obj,posgp)
            ngaus = size(posgp,2);
%             nelem = size(posgp,3);
            nelem = 2;
            x = posgp(1,:,:);
            y = posgp(2,:,:);
            fEl = obj.finiteElement;
            obj.deriv = zeros(2,length(fEl.shapeFunctions),ngaus,nelem);
            
            for s = 1:length(fEl.shapeFunctions)
                for i = 1:ngaus
                    obj.deriv(1,s,i,:) = fEl.shapeFunctionsDiff{s,1}(x(i),y(i));
                    obj.deriv(2,s,i,:) = fEl.shapeFunctionsDiff{s,2}(x(i),y(i));
                end
            end
        end
        
        function computeCases(obj)
            obj.iteration = [1 2 3;
                             2 3 1];
            obj.cases(:,:,1) = [1 4 5;
                4 2 3;
                5 4 3];
            obj.cases(:,:,2) = [1 4 3;
                4 2 5;
                4 5 3];
            obj.cases(:,:,3) = [1 4 5;
                1 2 4;
                5 4 3];
            obj.selectcases =  [1     0;
                2     0;
                3     3;
                0     2;
                0     1];
        end
        
    end
    
end
