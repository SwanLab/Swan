classdef ShFunc_Compliance < handle

    properties (Access = public)
        value
        gradient
    end

    properties (Access = private)
        Msmooth
        value0
    end

    properties (Access = private)
        mesh
        filter
        physicalProblem
        adjointProblem
        materialDerivative
    end

    methods (Access = public)

        function obj = ShFunc_Compliance(cParams)
            obj.init(cParams);
            obj.computeSmoothMassMatrix();
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
            obj.materialDerivative = cParams.materialDerivative;
        end

        function computeSmoothMassMatrix(obj)
            s.type  = 'MassMatrix';
            s.mesh  = obj.mesh;
            s.test  = P1Function.create(obj.mesh, 1);
            s.trial = P1Function.create(obj.mesh, 1);
            s.quadratureOrder = 'QUADRATICMASS';
            LHS = LHSintegrator.create(s);
            M   = LHS.compute();
            obj.Msmooth = M;
        end

        function computeFunction(obj)
            phy = obj.physicalProblem;
            stress = phy.stressFun.fValues;
            strain = phy.strainFun.fValues;
            dvolum = obj.mesh.computeDvolume(phy.strainFun.quadrature)';
            ngaus  = size(strain,2);
            nelem  = size(strain,3);

            c = zeros(nelem,ngaus);
            for igaus = 1:ngaus
                stressG = squeeze(stress(:,igaus,:));
                strainG  = squeeze(strain(:,igaus,:));
                e = stressG.*strainG;
                c(:,igaus) = c(:,igaus) + sum(e)';
            end
            int = c.*dvolum;
            obj.value = sum(int(:));
        end
        
        function computeGradient(obj)
            phy   = obj.physicalProblem;
            adj   = obj.adjointProblem;
            eu    = phy.strainFun.fValues;
            ep    = adj.strainFun.fValues;
            q     = phy.strainFun.quadrature;
            dC    = obj.materialDerivative.evaluate(q.posgp);
            nstre = size(eu,1);
            ngaus = size(eu,2);
            nelem = size(eu,3);
            g     = zeros(nelem,ngaus);
            for igaus = 1:ngaus
                for istre = 1:nstre
                    for jstre = 1:nstre
                        eui = squeeze(eu(istre,igaus,:));
                        epj = squeeze(ep(jstre,igaus,:));
                        dCij = squeeze(dC(istre,jstre,igaus,:));
                        g(:,igaus) = g(:,igaus) + (-eui.*dCij.*epj);
                    end
                end
            end
            obj.gradient = g;
        end

        function filterGradient(obj)
            g     = obj.gradient;
            nelem = size(g,1);
            ngaus = size(g,2);
            q     = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('LINEAR');
            s.fValues    = reshape(g',[1,ngaus,nelem]);
            s.mesh       = obj.mesh;
            s.quadrature = q;
            f            = FGaussDiscontinuousFunction(s);
            gradP1       = obj.filter.compute(f,'LINEAR');
            gf           = gradP1.fValues;
            gf           = obj.Msmooth*gf;
            s.fValues    = gf;
            obj.gradient = P1Function(s);
        end
    end
end