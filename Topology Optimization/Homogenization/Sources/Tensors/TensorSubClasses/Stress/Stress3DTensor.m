classdef Stress3DTensor < SecondOrder3DTensor ...
                            & StressDescriptor ...
    
    methods (Access = public)
        
        function obj = Stress3DTensor()
        end
        
        function makeItPlaneStressCompatible(obj)
            t = obj.tensorValue;
            psIndex = PlaneStressIndex();
            outPlane = psIndex.getOutPlaneIndex();
            conv     = TensorVoigtIndexTransformer3D();
            d = length(outPlane);
            for i = 1:d
                iv = outPlane(i);
                [it,jt] = conv.voigt2tensor(iv);
                t(it,jt) = 0;
                t(jt,it) = 0;
            end
            obj.tensorValue = t;
        end
    end
    
    
end

