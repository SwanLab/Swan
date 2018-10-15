classdef Filter_PDE_LevelSet_2D_Interior < Filter_PDE_LevelSet & Filter_LevelSet_2D_Interior
    methods
        function obj = Filter_PDE_LevelSet_2D_Interior(problemID,scale)
            obj@Filter_PDE_LevelSet(problemID,scale);
        end
        
         function preProcess(obj)
            preProcess@Filter_PDE(obj)
            preProcess@Filter_LevelSet_2D_Interior(obj)
        end
    end
end