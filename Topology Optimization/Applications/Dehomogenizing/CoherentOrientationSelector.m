classdef CoherentOrientationSelector < handle
    
    properties (Access = private)
       mesh
       orientation
    end
    
    methods (Access = public)
        
        function obj = CoherentOrientationSelector(cParams)
            obj.init(cParams);
        end
        
        function isCoherent = isOrientationCoherent(obj)
            nnode   = obj.mesh.nnodeElem;
            nElem   = obj.mesh.nelem;
            isCoh   = false(1,nnode,nElem);
            orient  = obj.orientation.fValues;            
            bN1 = squeeze(orient(:,1,:));
            for iNode = 1:nnode
                bNi = squeeze(orient(:,iNode,:));
                bN1bNI= dot(bN1,bNi);
                isCoh(1,iNode,:) = (bN1bNI)>0;
            end
            s.fValues  = isCoh;
            s.mesh     = obj.mesh;
            isCoherent = P1DiscontinuousFunction(s);            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh        = cParams.mesh;
            obj.orientation = cParams.orientation;
        end
        
    end
    
end