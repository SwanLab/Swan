classdef Filter_P1_LevelSet_2D_Interior < Filter_P1_LevelSet & Filter_LevelSet_2D_Interior
    methods        
         function preProcess(obj)
            preProcess@Filter_P1(obj)
            preProcess@Filter_LevelSet_2D_Interior(obj)
        end
    end
end