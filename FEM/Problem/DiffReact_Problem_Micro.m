classdef DiffReact_Problem_Micro < DiffReact_Problem

    methods (Access = protected)
        
        function setElement(obj)
            isRobinTermAdded = obj.isRobinTermAdded;
            bcType = obj.bcApplierType;            
            obj.element = Element_DiffReact_Micro(obj.mesh,obj.geometry,...
                obj.material,obj.dof,obj.problemData.scale,isRobinTermAdded,bcType,obj.interp);
        end
        
        function setDOFs(obj)
            obj.dof = DOF_DiffReact_Micro(obj.problemData.fileName,obj.geometry,obj.interp);
        end
        
        function setScale(obj)
            obj.problemData.scale = 'MICRO';
        end
        
    end
    
end
