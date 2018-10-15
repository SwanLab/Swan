classdef Filter_LevelSet_Interior < Filter_LevelSet
    
    properties (Access = private)
        shapeValues_FullCells
    end
    
    methods
        function preProcess(obj)
            preProcess@Filter_LevelSet(obj);
            obj.shapeValues_FullCells = obj.integrateFullCells;
        end
        
        function M2 = computeRHS(obj,x)
            obj.unfitted_mesh.computeMesh(x);
            obj.unfitted_mesh.computeDvoluCut;
            
            posgp_iso = obj.computePosGP(obj.unfitted_mesh.coord_iso_per_cell,obj.interpolation_unfitted,obj.quadrature_unfitted);
            obj.geometry.interpolation.computeShapeDeriv(squeeze(posgp_iso));
            
            shapeValues_CutCells = obj.integrateCutCells(obj.unfitted_mesh.cell_containing_subcell,obj.unfitted_mesh.dvolu_cut);
            shapeValues_All = obj.assembleShapeValues(obj.shapeValues_FullCells,shapeValues_CutCells);
            
            M2 = obj.rearrangeOutputRHS(shapeValues_All);
        end
        
        function S = computeVolume(obj,x)
            M2 = obj.computeRHS(x);
            S = sum(M2);
        end
    end
    
    methods (Access = private)
        function shapeValues_FullCells = integrateFullCells(obj)
            shapeValues_FullCells = zeros(size(obj.mesh.connec,1),size(obj.mesh.connec,2));
            for igauss = 1:size(obj.geometry.interpolation.shape,2)
                shapeValues_FullCells = shapeValues_FullCells + obj.geometry.interpolation.shape(:,igauss)'.*obj.geometry.dvolu(:,igauss);
            end
        end
        
        function shapeValued_CutCells = integrateCutCells(obj,containing_cell,dvolu_cut)
            dvolu_frac = sum(obj.geometry.dvolu,2)/obj.geometry.interpolation.dvolu;
            shapeValued_CutCells = obj.geometry.interpolation.shape'.*dvolu_cut.*dvolu_frac(containing_cell);
        end
        
        function shapeValues_AllCells = assembleShapeValues(obj,shapeValues_FullCells,shapeValues_CutCells)
            shapeValues_AllCells = zeros(size(obj.mesh.connec,1),size(obj.mesh.connec,2));
            shapeValues_AllCells(obj.unfitted_mesh.full_cells,:) = shapeValues_FullCells(obj.unfitted_mesh.full_cells,:);
            
            for i_subcell=1:size(shapeValues_CutCells,2)
                shapeValues_AllCells(:,i_subcell) = shapeValues_AllCells(:,i_subcell)+accumarray(obj.unfitted_mesh.cell_containing_subcell,shapeValues_CutCells(:,i_subcell),[obj.nelem,1],@sum,0);
            end
        end
    end
end