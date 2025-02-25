classdef testInverseOfInverseForStiffTensor < handle
    
    properties (Access = public)
        tol = 1e-10;
    end
    
    properties (Access = private)
        ctens
        invCtens
        invInvCtens
    end
    
    methods (Access = public)
        
        function obj = testInverseOfInverseForStiffTensor()
              obj.createStiffTensor()
              obj.createInvStiffTensor()
              obj.createInvInvStiffTensor()
        end
        
        function error = computeError(obj)
            c = obj.ctens.getValue();
            invInvC = obj.invInvCtens.getValue();
            error = norm(c(:) - invInvC(:));
        end

    end
    
    methods (Access = private)
        
        function createStiffTensor(obj)
            obj.ctens = Stiffness3DTensor();
            obj.ctens.createRandomTensor();
        end
        
        function createInvStiffTensor(obj)
            c = obj.ctens;
            obj.invCtens = Inverter.invert(c);
        end
        
        function createInvInvStiffTensor(obj)
            invC = obj.invCtens;
            obj.invInvCtens = Inverter.invert(invC);
        end

    end

end