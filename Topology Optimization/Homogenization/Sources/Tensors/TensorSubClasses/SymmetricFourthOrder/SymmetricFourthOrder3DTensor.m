classdef SymmetricFourthOrder3DTensor < FourthOrder3DTensor
    
    methods (Access = public)
        
        function obj = SymmetricFourthOrder3DTensor()
        end
        
        function createRandomTensor(obj)
            obj.createRandomTensor@FourthOrder3DTensor();
            obj.MakeMajorAndMinorSymmetrization();
        end
        
        function MakeMajorAndMinorSymmetrization(obj)
            Symmetrizer = FourthOrderSymmetrizer;
            obj.tensorValue = Symmetrizer.symmetrize(obj.tensorValue);
        end
        
    end
    
    
end

