classdef Filter_LevelSet_2D < Filter_LevelSet
    methods (Static)        
        function djacob = mapping(points,dvolu)
            v = diff(points);
            L = norm(v);
            djacob = L/dvolu;
        end
    end
end

