classdef Filter_PDE_LevelSet_3D_Boundary < Filter_PDE_LevelSet & Filter_LevelSet_3D_Boundary
    methods
        function preProcess(obj)
            preProcess@Filter_PDE(obj)
            preProcess@Filter_LevelSet_3D_Boundary(obj)
        end
    end
end