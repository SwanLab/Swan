classdef Integrator_Interior < Integrator  
    
    methods (Access = protected)
        
        function A = computeIntegral(obj,F1)
            if obj.isLeveSetCuttingMesh()
                shapeValues_CutCells = obj.integrateCutCells(F1);
                shapeValues_FullCells = obj.integrateFullCells(F1);
                shapeValues_All = obj.assembleShapeValues(shapeValues_CutCells,shapeValues_FullCells);
            else
                shapeValues_All = obj.integrateFullCells(F1);
            end
            
            A = obj.rearrangeOutputRHS(shapeValues_All);
        end
        
        function A = computeLHS(obj,globalConnec,npnod)
            interpolation = Interpolation.create(obj.meshBackground,'LINEAR');
            quadrature = obj.computeQuadrature(obj.meshBackground.geometryType);
            interpolation.computeShapeDeriv(quadrature.posgp);
            geometry = Geometry(obj.meshBackground,'LINEAR');
            geometry.computeGeometry(quadrature,interpolation);
            nelem = obj.meshBackground.nelem;
            Ae = zeros(interpolation.nnode,interpolation.nnode,nelem);            
            for igaus = 1:quadrature.ngaus
                for inode = 1:interpolation.nnode
                    for jnode = 1:interpolation.nnode
                        Ae(inode,jnode,:) = squeeze(Ae(inode,jnode,:)) + quadrature.weigp(igaus)*interpolation.shape(inode,igaus)...
                            *interpolation.shape(jnode,igaus)*geometry.djacob(:,igaus);
                    end
                end
            end
            
            nunkn1 = 1;
            nunkn2 = 1;
            nnode1 = size(globalConnec,2);
            nnode2 = size(globalConnec,2);
            idx1 = globalConnec';
            idx2 = globalConnec';
           
            A = sparse(npnod,npnod);
            for i = 1:nnode1*nunkn1
                for j = 1:nnode2*nunkn2
                    a = squeeze(Ae(i,j,:));
                    A = A + sparse(idx1(i,:),idx2(j,:),a,npnod,npnod);
                end
            end                
           
        end
        
    end
    
    methods (Access = private)
        
        function shapeValues_FullCells = integrateFullCells(obj,F1)            
            interpolation = Interpolation.create(obj.meshBackground,'LINEAR');
            quadrature = obj.computeQuadrature(obj.meshBackground.geometryType);
            interpolation.computeShapeDeriv(quadrature.posgp);
            geometry = Geometry(obj.meshBackground,'LINEAR');
            geometry.computeGeometry(quadrature,interpolation);
            
            shapeValues_FullCells = zeros(size(obj.meshBackground.connec));
            for igauss = 1:quadrature.ngaus
                shapeValues_FullCells = shapeValues_FullCells + interpolation.shape(:,igauss)'.*geometry.dvolu(:,igauss);
            end
        end
        
        function shapeValues_AllCells = assembleShapeValues(obj,shapeValues_CutCells,shapeValues_FullCells)
            interpolation = Interpolation.create(obj.meshBackground,'LINEAR');
            shapeValues_AllCells = zeros(size(obj.meshBackground.connec));
            shapeValues_AllCells(obj.meshUnfitted.backgroundFullCells,:) = shapeValues_FullCells(obj.meshUnfitted.backgroundFullCells,:);
            
            for i_subcell = 1:size(shapeValues_CutCells,2)
                shapeValues_AllCells(:,i_subcell) = shapeValues_AllCells(:,i_subcell)+accumarray(obj.meshUnfitted.cellContainingSubcell,shapeValues_CutCells(:,i_subcell),[interpolation.nelem,1],@sum,0);
            end
        end
        
    end
    
end

