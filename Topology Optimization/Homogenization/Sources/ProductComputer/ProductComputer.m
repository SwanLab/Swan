classdef ProductComputer < handle
    
    properties (Access = protected)
        secondOrder
        fourthOrder
        secondOrderOut
    end
    
    methods (Access = public, Static)
        
        function s = compute(C,e)
            factory = ProductComputerFactory();
            productComputer = factory.create(C,e);
            s = productComputer.getOutputSecondOrder();
        end
        
    end
    
    methods (Access = protected)
        
        function generate(obj,C,e)
            obj.init(C,e)
            obj.createOutSecondOrder()
            obj.computeProduct()
        end
        
    end
    
    
    methods (Access = private)
        
        function init(obj,C,e)
            obj.fourthOrder = C;
            obj.secondOrder = e;
        end
        
        
        function createOutSecondOrder(obj)
            C = obj.fourthOrder;
            e = obj.secondOrder;
            outCreator = ProductComputerOutputCreator(C,e);
            obj.secondOrderOut = outCreator.getSecondOrderOut;
        end
        
        function s = getOutputSecondOrder(obj)
            s = obj.secondOrderOut;
        end
    end
    
    
    methods (Access = protected, Abstract)
        computeProduct(obj)
    end
    
end

