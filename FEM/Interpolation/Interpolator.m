classdef Interpolator < handle
    
    properties (Access = private)
       sMesh
       interpolation
       cellFinder
       zGrid
       zInterp
       zInterpDeriv       
    end    
    
    methods (Access = public)
    
        function obj = Interpolator(cParams)
            obj.init(cParams);
            obj.createInterpolation();
        end
        
        function setValues(obj,x,y)
            s.mesh     = obj.sMesh;
            s.points.x = x;
            s.points.y = y;
            obj.cellFinder = CellFinderInStructuredMesh(s);            
            obj.evaluateShapeFunctions();
        end
        
        function [z,dz] = interpolate(obj,z)
            obj.zGrid  = z;                 
            obj.obtainInterpolationValues();
            obj.obtainInterpolationDerivativesValues();
            z = obj.zInterp;
            dz = obj.zInterpDeriv;
        end
        
    end
    
    
    methods (Access = private)
                
        function init(obj,cParams)
            obj.sMesh = cParams.mesh;                    
        end
        
        function createInterpolation(obj)
            m = obj.sMesh.mesh;
            int = Interpolation.create(m,'LINEAR');
            obj.interpolation = int;    
        end        
        
        function evaluateShapeFunctions(obj)
            nC = obj.cellFinder.naturalCoord;
            obj.interpolation.computeShapeDeriv(nC);
        end
                
        function obtainInterpolationValues(obj)
            shapes = obj.interpolation.shape;        
            z(:,1) = reshape(obj.zGrid',[],1);
            v      = zeros(size(shapes,2),1);
            for inode = 1:size(shapes,1)
                nodes = obj.cellFinder.cells(:,inode);
                znode = z(nodes);
                v(:,1) = v(:,1) + znode.*shapes(inode,:)';
            end            
            obj.zInterp = v;
        end
        
        function obtainInterpolationDerivativesValues(obj)
            shapeDeriv = obj.interpolation.deriv;            
            z(:,1) = reshape(obj.zGrid',[],1);
            v      = zeros(size(shapeDeriv,3),size(shapeDeriv,1));
            for inode = 1:size(shapeDeriv,2)
                nodes = obj.cellFinder.cells(:,inode);
                znode = z(nodes);
                sh = squeeze(shapeDeriv(:,inode,:));
                v = v + bsxfun(@(x,y) x.*y,znode,sh');
            end            
            obj.zInterpDeriv = v;
        end
        

    end
    
    
end