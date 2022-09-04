classdef Projector_P0toP1 < handle
    
    properties (Access = private)
        
    end
    
    methods (Access = public)

        function obj = Projector_P0toP1(cParams)
            obj.init(cParams);
        end

        function project(obj)
            obj.computeLHS();
            obj.computeRHS();
            obj.solve();
        end

    end

    methods (Access = private)
        
        function init(obj, cParams)
        end

        function computeLHS(obj)
        end

        function computeRHS(obj)
            % ez
        end

        function solve(obj)
        end
        
    end

end

