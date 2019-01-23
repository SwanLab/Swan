classdef Filter_PDE_LevelSet_3D_Interior < Filter_PDE_LevelSet & Filter_LevelSet_3D_Interior
    methods       
         function preProcess(obj)
            preProcess@Filter_PDE(obj)
            preProcess@Filter_LevelSet_3D_Interior(obj)
        end
    end
end