classdef testPlotting_Composite < testPlotting
    methods (Access = protected)
        function hasPassed = hasPassed(obj)
            coordsOk = obj.checkCoordinates(obj.storedVar{1});
            connecsOk = obj.checkConnectivities(obj.storedVar{2});
            hasPassed = coordsOk && connecsOk;
        end
    end
    
    methods (Access = private)
        function coordsOk = checkCoordinates(obj,storedCoords)
            boxMeshesCoords = cell(1,obj.mesh.nActiveBoxFaces);
            for iactive = 1:obj.mesh.nActiveBoxFaces
                iface = obj.mesh.activeBoxFaceMeshesList(iactive);
                boxMeshesCoords{iactive} = obj.mesh.boxFaceMeshes{iface}.coord;
            end
            coordsOk = obj.checkComposedVar(obj.mesh.innerMeshOLD.coord,boxMeshesCoords,storedCoords);
        end
        
        function connecsOk = checkConnectivities(obj,storedConnecs)
            boxMeshesConnecs = cell(1,obj.mesh.nActiveBoxFaces);
            for iactive = 1:obj.mesh.nActiveBoxFaces
                iface = obj.mesh.activeBoxFaceMeshesList(iactive);
                boxMeshesConnecs{iactive} = obj.mesh.boxFaceMeshes{iface}.connec;
            end
            connecsOk = obj.checkComposedVar(obj.mesh.innerMeshOLD.connec,boxMeshesConnecs,storedConnecs);
        end
    end
    
    methods (Access = private)
        function varOk = checkComposedVar(obj,interiorVar,boxMeshesVar,storedVar)
            interiorVarOk = obj.checkVar(interiorVar,storedVar{1});
            boxMeshesVarOk = true;
            for iactive = 1:obj.mesh.nActiveBoxFaces
                boxMeshesVarOk = boxMeshesVarOk && obj.checkVar(boxMeshesVar{iactive},storedVar{1+iactive});
            end
            varOk = interiorVarOk && boxMeshesVarOk;
        end
    end
    
    methods (Access = private, Static)
        function varOk = checkVar(var,ref)
            if all(size(var) == size(ref))
                if all(all(var == ref))
                    varOk = true;
                else
                    varOk = false;
                end
            else
                varOk = false;
            end
        end
    end
end

