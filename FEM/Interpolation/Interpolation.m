classdef Interpolation < handle
    
    properties (GetAccess = public, SetAccess = protected)
      

        order
        
        ndime
        nnode       
        
        pos_nodes
        shape
        deriv
        isoDv
        
        iteration
        cases
        selectcases
        main_loop
        extra_cases
    end
    
    properties (Access = protected)
        type
    end
    
    methods (Static, Access = public)
        
        function obj = create(mesh,order)
            cParams.mesh = mesh;
            cParams.order = order;
            f = InterpolationFactory;
            obj = f.create(cParams);
        end
        
    end
    
    methods (Access = public) 
        
        function fxV = interpolateFunction(obj,xV,func)
            obj.computeShapeDeriv(xV);
            shapes = obj.shape;
            nNode  = size(shapes,1);
            nGaus  = size(shapes,2); 
            nF     = size(func,1);                        
            nElem  = size(func,3);
            fxV = zeros(nF,nGaus,nElem);
            for kNode = 1:nNode
                shapeKJ = shapes(kNode,:,:);
                fKJ     = func(:,kNode,:);
                f = bsxfun(@times,shapeKJ,fKJ);
                fxV = fxV + f;
            end
        end        
        
        
    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
            obj.type  = cParams.mesh.geometryType;
            obj.order = cParams.order;            
        end
        
    end
    
    methods (Abstract)
        computeShapeDeriv(obj)
    end
end
