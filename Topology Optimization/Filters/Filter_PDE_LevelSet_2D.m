classdef Filter_PDE_LevelSet_2D < Filter_PDE_LevelSet & Filter_LevelSet_2D
    methods
        function obj = Filter_PDE_LevelSet_2D(problemID,scale)
            obj@Filter_PDE_LevelSet(problemID,scale);
        end
        
        function preProcess(obj)
            preProcess@Filter_PDE(obj)
            preProcess@Filter_LevelSet_2D(obj)
        end
    end
end