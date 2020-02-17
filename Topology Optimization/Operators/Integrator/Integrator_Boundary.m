classdef Integrator_Boundary < IntegratorUnfitted
    
    properties (Access = private)
        shapes
        cutShapes
    end
    
    methods (Access = public)
        
        function A = integrate(obj,F1)
            obj.initShapes();
            % if obj.isLeveSetCuttingMesh()
            obj.cutShapes = obj.evaluateCutShapes(F1);
            obj.assembleShapes();
            % else
            %    a = 1
            % end
            A = obj.rearrangeOutputRHS(obj.shapes);
        end
        
        
    end
    
    methods (Access = private)
        
        function initShapes(obj)
            nelem = obj.meshBackground.nelem;
            nnode = obj.meshBackground.nnode;
            obj.shapes = zeros(nelem,nnode);
        end
        
        function assembleShapes(obj)
            obj.assembleCutShapes();
        end
        
        function assembleCutShapes(obj)
            nelem = obj.meshBackground.nelem;
            cell = obj.meshUnfitted.cellContainingSubcell;
            nnode = obj.meshBackground.nnode;
            
            for iNode = 1:nnode
                csNode = obj.cutShapes(:,iNode);
                csNodeGlobal  = accumarray(cell,csNode,[nelem,1],@sum,0);
                sNode  = obj.shapes(:,iNode);
                obj.shapes(:,iNode) = sNode + csNodeGlobal;
            end
        end
        
    end
    
end

