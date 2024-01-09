classdef ShFunc_Compliance < handle

    properties (Access = public)
        value
        gradient
    end

    properties (Access = private)
        mesh
        filter
        physicalProblem
        adjointProblem
        material
        materialDerivative
    end

    methods (Access = public)

        function obj = ShFunc_Compliance(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            obj.computeFunction();
            obj.computeGradient();
            obj.filterGradient();
        end
    end
    
    methods (Access = private)

        function init(obj,cParams)
            obj.mesh               = cParams.mesh;
            obj.filter             = cParams.filter;
            obj.physicalProblem    = cParams.physicalProblem;
            obj.adjointProblem     = cParams.physicalProblem;
            obj.material           = cParams.material;
            obj.materialDerivative = cParams.materialDerivative;
        end

        function computeFunction(obj)
            u      = obj.physicalProblem.uFun;
            q      = obj.physicalProblem.getQuadrature();
            Cij    = obj.material.evaluate(q.posgp);
            strain = u.computeSymmetricGradient(q);
            strain.applyVoigtNotation();
            strainj(:,1,:,:) = strain.fValues;
            stressj          = pagemtimes(Cij,strainj);
            stressj          = permute(stressj, [1 3 4 2]);
            s.quadrature     = q;
            s.fValues        = stressj;
            s.mesh           = obj.mesh;
            stress           = FGaussDiscontinuousFunction(s);
            iPar.mesh        = obj.mesh;
            iPar.quadType    = q.order;
            int              = IntegratorScalarProduct(iPar);
            obj.value        = int.compute(strain,stress);
        end

        function computeGradient(obj)
            q            = obj.physicalProblem.getQuadrature();
            u            = obj.physicalProblem.uFun;
            p            = obj.adjointProblem.uFun;
            eu           = u.computeSymmetricGradient(q);
            eu.applyVoigtNotation();
            ep           = p.computeSymmetricGradient(q);
            ep.applyVoigtNotation();
            dC           = obj.materialDerivative.evaluate(q.posgp);
            epi(1,:,:,:) = ep.fValues;
            euj(:,1,:,:) = eu.fValues;
            dCeuj        = pagemtimes(dC,euj);
            contGrad     = -pagemtimes(epi,dCeuj);
            contGrad     = squeezeParticular(contGrad,1);
            s.fValues    = contGrad;
            s.mesh       = obj.mesh;
            s.quadrature = q;
            g            = FGaussDiscontinuousFunction(s);
            ss.mesh      = obj.mesh;
            ss.type      = 'ShapeFunction';
            ss.quadType  = q.order;
            int          = RHSintegrator.create(ss);
            test         = P1Function.create(obj.mesh,1);
            Dxc          = int.compute(g,test);
            s.fValues    = Dxc;
            obj.gradient = P1Function(s);
        end

        function filterGradient(obj)
            g       = obj.gradient;
            regGrad = obj.filter.compute(g,'LINEAR');
            obj.gradient = regGrad;
        end
    end
end