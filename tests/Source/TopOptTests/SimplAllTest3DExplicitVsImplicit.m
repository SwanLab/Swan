classdef SimplAllTest3DExplicitVsImplicit < SimplAllTestExplicitVsImplicit
     
    methods (Access = public)
       
        function obj = SimplAllTest3DExplicitVsImplicit()
            obj.dim = '3D';
            obj.init();
            obj.createProperties();
            obj.computeExplicitMaterialInterpolation();
            obj.computeImplicitMaterialInterpolation();            
        end
        
    end

end