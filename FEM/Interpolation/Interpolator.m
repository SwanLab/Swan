classdef Interpolator < handle
    
    properties (Access = private)
       sMesh
       interpolation
       cellFinder
       zGrid
       zValues
       nComponents
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
            obj.nComponents = size(z,3);
            obj.zValues  = obj.computeZvalues(z);                 
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
            v = obj.interpolateValues(shapes,obj.zValues);
            obj.zInterp = v;
        end
        
        function v = interpolateValues(obj,shapes,f)
            shapesT = shapes';            
            v = bsxfun(@times,shapesT,f);
            v = permute(v,[3 1 2]);
            v = sum(v,3);
        end
        
        function zT = computeZvalues(obj,z)
            z = permute(z,[2 1 3]);
            z = reshape(z,[],obj.nComponents);
            nodes = obj.cellFinder.cells;
            [nelem,nnode] = size(nodes);
            znode = z(nodes(:),:);
            zT = reshape(znode,nelem,nnode,obj.nComponents);    
        end
                
        function obtainInterpolationDerivativesValues(obj)
            shapeDeriv = obj.interpolation.deriv;            
            ndir  = size(shapeDeriv,1);
            nelem = size(shapeDeriv,3);
            v = zeros(obj.nComponents,nelem,ndir);            
            for idir = 1:ndir
                shapes = squeeze(shapeDeriv(idir,:,:));
                vi = obj.interpolateValues(shapes,obj.zValues);
                v(:,:,idir) = vi;
            end
            obj.zInterpDeriv = v;
        end         
        
    end
    
    
end