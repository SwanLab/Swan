classdef SimplAllTest2DExplicitVsImplicit < SimplAllTestExplicitVsImplicit
     
    methods (Access = public)
       
        function obj = SimplAllTest2DExplicitVsImplicit()
            obj.dim = '2D';
            obj.init();
            obj.createProperties();
            obj.computeExplicitMaterialInterpolation();
            obj.computeImplicitMaterialInterpolation();            
        end
        
    end

end