classdef Filter_P1_LevelSet_2D_Interior < Filter_P1_LevelSet & Filter_LevelSet_2D_Interior
    methods
        function obj = Filter_P1_LevelSet_2D_Interior(problemID,scale)
            obj@Filter_P1_LevelSet(problemID,scale);
        end
        
         function preProcess(obj)
            preProcess@Filter_P1(obj)
            preProcess@Filter_LevelSet_2D_Interior(obj)
        end
    end
end