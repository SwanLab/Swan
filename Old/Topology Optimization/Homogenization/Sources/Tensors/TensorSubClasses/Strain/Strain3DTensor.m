classdef Strain3DTensor < SecondOrder3DTensor  ...
                          & StrainDescriptor
    
    methods (Access = public)
        
        function obj = Strain3DTensor()
        end
        
        function makeItPlaneStressCompatible(obj,Ch)
            stress = Stress3DTensor();
            stress.createRandomTensor();
            stress.makeItPlaneStressCompatible();
            Chinv = Inverter.invert(Ch);
            strain = ProductComputer.compute(Chinv,stress);
            obj.setValue(strain.getValue())
        end
    end

end

