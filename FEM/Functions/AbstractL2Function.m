classdef AbstractL2Function < L2Function
       
   properties (Access = public)
        ndimf
    end
    
    properties (Access = private)
        evaluateF
    end
    
    methods (Access = public)

        function obj = AbstractL2Function(cParams)
            obj.mesh      = cParams.mesh;
            obj.evaluateF = cParams.evaluate;
            obj.ndimf     = cParams.ndimf;
        end

        function f = evaluate(obj,xV)
           f = obj.evaluateF(xV);
        end
    end

end