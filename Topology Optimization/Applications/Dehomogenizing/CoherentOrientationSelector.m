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
%             connec  = obj.mesh.connec;
            orient  = obj.orientation.fValues;
            isCoh   = false(nElem,nnode);
%             for iElem = 1:nElem
%                 nodeRef = connec(iElem,1);
%                 o1 = orient(nodeRef,:);
%                 for iNode = 1:nnode
%                     nodeI = connec(iElem,iNode);
%                     oI   = orient(nodeI,:);
%                     o1oI = dot(o1,oI);
%                     isCoh(iElem,iNode) = (o1oI)>0;
%                 end
%             end
            for iElem = 1:nElem % !!
                o1 = orient(:,1,iElem);
                for iNode = 1:nnode
                    oI = orient(:,iNode,iElem);
                    o1oI = dot(o1,oI);
                    isCoh(iElem,iNode) = (o1oI)>0;
                end
            end
            
            a.fValues = permute(isCoh, [3 2 1]);
            a.mesh = obj.mesh;
            isCoherent = P1DiscontinuousFunction(a);            
           % isCoherent = isCoh;

        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh        = cParams.mesh;
            obj.orientation = cParams.orientation;
        end
        
    end
    
end