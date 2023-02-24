classdef OrientationVectors < handle
    
    properties (GetAccess = public, SetAccess = private)
        value
        isCoherent
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
      mesh
      theta
    end
    
    methods (Access = public)
        
        function obj = OrientationVectors(cParams)
            obj.init(cParams)
            obj.createOrientationVector();
            obj.computeIsOrientationCoherent();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh  = cParams.mesh;
            obj.theta = cParams.theta;
        end
        
        function createOrientationVector(obj)
            a1(:,1) = cos(obj.theta);
            a1(:,2) = sin(obj.theta);
            a2(:,1) = -sin(obj.theta);
            a2(:,2) = cos(obj.theta);
            a(:,:,1) = a1;
            a(:,:,2) = a2;
            nDim = obj.mesh.ndim;
            orientation = cell(nDim,1);
            for iDim = 1:nDim
                s.fValues = a(:,:,iDim);
                s.mesh   = obj.mesh;
                bf = P1Function(s);
                orientation{iDim} = bf;
            end 
            obj.value = orientation;
        end
        
        function computeIsOrientationCoherent(obj)
            s.mesh        = obj.mesh;
            s.orientation = obj.value{1}.project('P1D');
            c = CoherentOrientationSelector(s);
            isC = c.isOrientationCoherent();
            obj.isCoherent = isC;            
        end
        
    end
    
end