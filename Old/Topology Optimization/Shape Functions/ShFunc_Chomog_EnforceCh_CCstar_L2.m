classdef ShFunc_Chomog_EnforceCh_CCstar_L2 < ShFunc_Chomog_EnforceCh
   
    methods (Access = public)
        
        function obj=ShFunc_Chomog_EnforceCh_CCstar_L2(cParams)
            obj.initChomog(cParams);
            obj.computeChTarget(cParams.ChTarget);
            obj.computeC0();
            obj.computeWeights();
            obj.pNorm = 2;
        end
        
    end
end
