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
            obj.ctens = FourthOrderTensor();
            obj.ctens.createRandomTensor();            
        end
        
        function createInvStiffTensor(obj)
            c = obj.ctens;
            obj.invCtens = obj.invertTensor(c);
        end
        
        function createInvInvStiffTensor(obj)            
            invC = obj.invCtens;
            obj.invInvCtens = obj.invertTensor(invC);
        end

    end
    
    methods (Access = private, Static)

        function invCtensor = invertTensor(c)
            invC = Inverter.invert(c);            
            invCtensor = FourthOrderTensor();
            invCtensor.setValue(invC);
        end        
    end
    
    methods (Access = protected)
        function hasPassed = hasPassed(obj)
            c = obj.ctens.getValue();
            invInvC = obj.invInvCtens.getValue();
            hasPassed = norm(c(:) - invInvC(:)) < 1e-12;
        end
    end
end

