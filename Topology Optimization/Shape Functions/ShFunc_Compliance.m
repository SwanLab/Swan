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
            Cij    = squeeze(obj.material.evaluate(q.posgp));
            strain = u.computeSymmetricGradient(q);
            strain.applyVoigtNotation();
            strainj       = strain.fValues;
            stressj       = pagemtimes(Cij,strainj);
            s.quadrature  = q;
            s.fValues     = stressj;
            s.mesh        = obj.mesh;
            stress        = FGaussDiscontinuousFunction(s);
            iPar.mesh     = obj.mesh;
            iPar.quadType = q.order;
            int           = IntegratorScalarProduct(iPar);
            obj.value     = int.compute(strain,stress);
        end

        function computeGradient(obj)
            q            = obj.physicalProblem.getQuadrature();
            u            = obj.physicalProblem.uFun;
            p            = obj.adjointProblem.uFun;
            eu           = u.computeSymmetricGradient(q);
            eu.applyVoigtNotation();
            ep           = p.computeSymmetricGradient(q);
            ep.applyVoigtNotation();
            dC           = squeeze(obj.materialDerivative.evaluate(q.posgp));
            epi          = permute(ep.fValues,[2 1 3]);
            euj          = eu.fValues;
            dCeuj        = pagemtimes(dC,euj);
            contGrad     = -pagemtimes(epi,dCeuj);
            s.fValues    = contGrad;
            s.mesh       = obj.mesh;
            s.quadrature = q;
            g            = FGaussDiscontinuousFunction(s);
            ss.mesh      = obj.mesh;
            ss.type      = 'ShapeFunction';
            ss.quadType  = 'LINEAR';
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