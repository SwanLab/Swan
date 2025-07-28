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


        function lhs = computeElementalLHS(obj,f)
            nElem  = obj.mesh.nelem;
            lhs    = zeros(size(f,1),size(f,2),nElem);

            J = Jacobian(obj.mesh);
            %J    = obj.mesh.getJacobian();
            detJ = Det(J);
            %detJ = DetJ(obj.mesh);

            xV = obj.quadrature.posgp;
            w  = obj.quadrature.weigp;
            for i = 1:size(f,1)
                for j = 1:size(f,2)
                    int = (f{i,j}.*detJ)*w';
                    lhs(i,j,:) = lhs(i,j,:) + int.evaluate(xV);
                end
            end
        end

        function init(obj,cParams)
            obj.mesh  = cParams.mesh;
            obj.setQuadratureOrder(cParams);
        end

    end
end