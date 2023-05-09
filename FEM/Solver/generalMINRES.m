classdef generalMINRES < handle

    properties (Access = public)
        iter
        xPrevIt
        Yk
    end

    methods (Access = public)
        function obj = generalMINRES()
            obj.init();
        end

        function u = solve(A,b)

            if (obj.iter <= 5)
                solver = MINRES_Pol
                u = MINRES_Pol.solve(A,b);
            else 
                u = actualMINRES
        end
    end

    methods (Access = private)
        function init(obj)
            obj.iter = 1;
            obj.xPrevIt = zeros(20200,1); 
        end
    end
end