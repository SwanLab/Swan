classdef Filter_PDE_LevelSet_3D_Boundary < Filter_PDE_LevelSet & Filter_LevelSet_3D_Boundary
    methods
        function obj = Filter_PDE_LevelSet_3D_Boundary(problemID,scale)
            obj@Filter_PDE_LevelSet(problemID,scale);
        end
        
         function preProcess(obj)
            preProcess@Filter_PDE(obj)
            preProcess@Filter_LevelSet_3D_Boundary(obj)
        end
    end
end