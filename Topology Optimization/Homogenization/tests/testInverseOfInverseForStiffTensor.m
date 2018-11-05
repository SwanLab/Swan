classdef testInverseOfInverseForStiffTensor < test
    
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
        
    end
    
    methods (Access = private)
        
        function createStiffTensor(obj)
            obj.ctens = SymmetricFourthOrder3DTensor();
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
       
    methods (Access = protected)
        function hasPassed = hasPassed(obj)
            c = obj.ctens.getValue();
            invInvC = obj.invInvCtens.getValue();
            hasPassed = norm(c(:) - invInvC(:)) < 1e-10;
        end
    end
end

