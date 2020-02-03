classdef Integrator_Interior < IntegratorUnfitted
    
    properties (Access = private)
        shapes
        cutShapes
        innerShapes
    end
    
    methods (Access = public)
        
        function A = integrate(obj,F1)
            obj.initShapes();
            if obj.isLeveSetCuttingMesh()
                obj.cutShapes = obj.evaluateCutShapes(F1);
                obj.evaluateInnerShapes(F1);
                obj.assembleShapes();
            else
                obj.evaluateInnerShapes(F1);
                obj.assembleInnerShapes();
            end
            A = obj.rearrangeOutputRHS(obj.shapes);
        end
        
    end
    
    methods (Access = private)
        
        function evaluateInnerShapes(obj,F1)
            interpolation = Interpolation.create(obj.meshBackground,'LINEAR');
            quadrature = obj.computeQuadrature(obj.meshBackground.geometryType);
            interpolation.computeShapeDeriv(quadrature.posgp);
            geometry = Geometry(obj.meshBackground);
            geometry.computeGeometry(quadrature,interpolation);
            
            obj.innerShapes = zeros(size(obj.meshBackground.connec));
            for igauss = 1:quadrature.ngaus
                obj.innerShapes = obj.innerShapes + interpolation.shape(:,igauss)'.*geometry.dvolu(:,igauss);
            end
        end
        
        function initShapes(obj)
            nelem = obj.meshBackground.nelem;
            nnode = obj.meshBackground.nnode;
            obj.shapes = zeros(nelem,nnode);
        end
        
        function assembleShapes(obj)
            obj.assembleInnerShapes();
            obj.assembleCutShapes();
        end
        
        function assembleInnerShapes(obj)
            innerCells = obj.meshUnfitted.backgroundFullCells;
            obj.shapes(innerCells,:) = obj.innerShapes(innerCells,:);
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

