classdef LHSIntegratorStiffnessElastic < LHSIntegrator

    methods (Access = public)
        function obj = LHSIntegratorStiffnessElastic(cParams)
               obj@LHSIntegrator(cParams)
         %   obj.init(cParams);
        end



        function LHS = compute(obj,f,test,trial)
            lhs = obj.computeElementalLHS(f);
            LHS = obj.assembleMatrix(lhs,test,trial);
        end
    end

    methods (Access = protected)

        function init(obj,cParams)
            obj.mesh  = cParams.mesh;
            obj.setQuadratureOrder(cParams);
        end

    end
end