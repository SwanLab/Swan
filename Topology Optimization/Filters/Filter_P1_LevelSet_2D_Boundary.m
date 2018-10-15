classdef Filter_P1_LevelSet_2D_Boundary < Filter_P1_LevelSet & Filter_LevelSet_2D_Boundary
    methods
        function obj = Filter_P1_LevelSet_2D_Boundary(problemID,scale)
            obj@Filter_P1_LevelSet(problemID,scale);
        end
        
         function preProcess(obj)
            preProcess@Filter_P1(obj)
            preProcess@Filter_LevelSet_2D_Boundary(obj)
        end
    end
end