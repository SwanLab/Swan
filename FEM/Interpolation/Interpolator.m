classdef Interpolator < handle
    
    properties (Access = private)
       mesh
       interpolation
       cellFinder
       zGrid
       zInterp
       zInterpDer       
    end    
    
    methods (Access = public)
    
        function obj = Interpolator(cParams)
            obj.init(cParams);
            obj.createInterpolation();
        end
        
        function setValues(obj,x,y)
            s.mesh     = obj.mesh;
            s.points.x = x;
            s.points.y = y;
            obj.cellFinder = CellFinderInStructuredMesh(s);            
            obj.evaluateShapeFunctions();
        end
        
        function z = interpolate(obj,z)
            obj.zGrid  = z;                 
            obj.obtainInterpolationValues();
            %obj.obtainInterpolationDerivativesValues();
            z = obj.zInterp;
        end
        
        function z = interpolateDerivative(obj,z)
            
        end
        
    end
    
    
    methods (Access = private)
                
        function init(obj,cParams)
            obj.mesh = cParams.mesh;                    
        end
        
        function createInterpolation(obj)
            int = Interpolation.create(obj.mesh,'LINEAR');
            obj.interpolation = int;    
        end        
        
        function evaluateShapeFunctions(obj)
            nC = obj.cellFinder.naturalCoord;
            obj.interpolation.computeShapeDeriv(nC);
        end
                
        function obtainInterpolationValues(obj)
            shapes = obj.interpolation.shape;            
            obj.zInterp = obj.sumNodesContribution(shapes);
        end
        
        function obtainInterpolationDerivativesValues(obj)
            shapeDeriv = obj.interpolation.deriv;            
            obj.zInterpDer = obj.sumNodesContribution(shapeDeriv);
        end
        
        function v = sumNodesContribution(obj,shapes)
            z(:,1) = reshape(obj.zGrid',[],1);
            v      = zeros(size(shapes,2),1);
            for inode = 1:size(shapes,1)
                nodes = obj.cellFinder.cells(:,inode);
                znode = z(nodes);
                v(:,1) = v(:,1) + znode.*shapes(inode,:)';
            end
        end
        
    end
    
    
end