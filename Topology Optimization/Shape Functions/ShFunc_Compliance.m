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
            q     = obj.physicalProblem.getQuadrature();
            u     = obj.physicalProblem.uFun;
            p     = obj.adjointProblem.uFun;
            eu    = u.computeSymmetricGradient(q);
            eu.applyVoigtNotation();
            ep    = p.computeSymmetricGradient(q);
            ep.applyVoigtNotation();
            dC    = obj.materialDerivative.evaluate(q.posgp);
            nstre = size(eu.fValues,1);
            ngaus = size(eu.fValues,2);
            nelem = size(eu.fValues,3);
            g     = zeros(nelem,ngaus);
            for igaus = 1:ngaus
                for istre = 1:nstre
                    for jstre = 1:nstre
                        eui = squeeze(eu.fValues(istre,igaus,:));
                        epj = squeeze(ep.fValues(jstre,igaus,:));
                        dCij = squeeze(dC(istre,jstre,igaus,:));
                        g(:,igaus) = g(:,igaus) + (-eui.*dCij.*epj);
                    end
                end
            end
            s.fValues    = reshape(g',[1,ngaus,nelem]);
            s.mesh       = obj.mesh;
            s.quadrature = q;
            obj.gradient = FGaussDiscontinuousFunction(s);
        end

        function filterGradient(obj)
            g       = obj.gradient;
            regGrad = obj.filter.compute(g,'LINEAR');
            obj.gradient = regGrad;
        end
    end
end