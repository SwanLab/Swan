classdef ShFunc_Compliance < handle

    properties (Access = public)
        value
        gradient
    end

    properties (Access = private)
        mesh
        filter
        stateProblem
        adjointProblem
  %      material
   %     materialDerivative
    end

    methods (Access = public)

        function obj = ShFunc_Compliance(cParams) %change name
            obj.init(cParams);
        end

        function compute(obj,x)
            obj.computeFunction(x);
            obj.computeGradient(x);
        end
    end
    
    methods (Access = private)

        function init(obj,cParams)
            obj.mesh               = cParams.mesh;
            obj.filter             = cParams.filter;
            obj.stateProblem       = cParams.stateProblem; %state
            obj.adjointProblem     = cParams.stateProblem; 
            %obj.material           = cParams.material;
            %obj.materialDerivative = cParams.materialDerivative;
            %material Interpolator
        end

        function computeMaterial(obj,x)
            C = obj.createMaterial(x);
            obj.stateProblem.updateMaterial(C);
            obj.stateProblem.solve();
            Cij    = obj.material.evaluate(q.posgp);
        end

        function computeFunction(obj,x)
            C    = obj.computeMaterial(x);
            u      = obj.stateProblem.uFun;
            q      = obj.stateProblem.getQuadrature();% Constriur quad, elimate getQuadr

            % ComputeStrain
            strain = u.computeSymmetricGradient(q);
            strain.applyVoigtNotation(); %strain = strain.obtainVoigtFormat();
            strainj(:,1,:,:) = strain.fValues;
            % ComputeStress
            stressV          = pagemtimes(Cij,strainj);
            stressV          = permute(stressV, [1 3 4 2]);
            s.quadrature     = q;
            s.fValues        = stressV;
            s.mesh           = obj.mesh;
            stress           = FGaussDiscontinuousFunction(s); %Use create

            iPar.mesh        = obj.mesh; % Use s
            iPar.quadType    = q.order;
            int              = IntegratorScalarProduct(iPar);
            J                = int.compute(strain,stress);
            obj.value        = J;
        end

        function computeGradient(obj)
            dj = obj.computeDJ();
            g  = obj.filter.compute(dj,'LINEAR'); 
            obj.gradient = g;
        end

        function dj = computeDJ(obj)
            q            = obj.stateProblem.getQuadrature(); %Not any more
            u            = obj.stateProblem.uFun;
            p            = obj.adjointProblem.uFun;
            dCij           = obj.materialDerivative.evaluate(q.posgp);            
            %ComputeStateStrain
            eu           = u.computeSymmetricGradient(q);
            eu.applyVoigtNotation();
            %ComputeAdjointStrain
            ep           = p.computeSymmetricGradient(q);
            ep.applyVoigtNotation();
            epi(1,:,:,:) = ep.fValues;
            euj(:,1,:,:) = eu.fValues;
            dStress = pagemtimes(dCij,euj);
            dj      = pagemtimes(epi,dStress);
            dj      = squeezeParticular(-dj,1);                        
            s.fValues    = dj;
            s.mesh       = obj.mesh;
            s.quadrature = q;
            dj           = FGaussDiscontinuousFunction(s); % create 
        end
  
    end
end