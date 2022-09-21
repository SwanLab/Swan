classdef ProjectorToP1discont < handle
    % OBJECTIU: encara no fa falta resoldre el sistema matricial. Simplement
    % fer switch depenent de la funcio d'entrada (P0 o P1) i trobar "a lo
    % cutre" la P1discont. Ja mes endavant fem el cas general. Ja hi ha
    % fragments de codi escampats que resolen P0toP1discont i P1toP1discont.

    properties (Access = public)

    end

    properties (Access = private)
        mesh
        connec
        quadOrder
    end

    properties (Access = private)
        quadrature
        field
        meshD
    end

    methods (Access = public)

        function obj = ProjectorToP1discont(cParams)
            obj.init(cParams);
            obj.createDiscontinuousMesh();
            obj.createQuadrature();
            obj.createField();
        end

        function xFun = project(obj, x)
            LHS = obj.computeLHS();
%             RHS = obj.computeRHS(x);
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.mesh   = cParams.mesh;
            obj.connec = cParams.connec;
            obj.quadOrder = 'LINEAR';
        end

        function createDiscontinuousMesh(obj)
            obj.meshD = obj.mesh.createDiscontinuousMesh();
        end

        function q = createQuadrature(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature(obj.quadOrder);
            obj.quadrature = q;
        end

        function createField(obj)
            s.mesh               = obj.meshD;
            s.ndimf              = 1;
            s.interpolationOrder = 'LINEAR';
            s.quadratureOrder    = 'QUADRATIC';
            obj.field = Field(s);
        end
        
        function LHS = computeLHS(obj)
            s.type  = 'MassMatrix';
            s.mesh  = obj.meshD;
            s.field = obj.field;
            lhs = LHSintegrator.create(s);
            LHS = lhs.compute();
        end

        function RHS = computeRHS(obj, fun)
            fDisc = fun.computeDiscontinuousField();
            fVals = fDisc.fValues;
            dV = obj.mesh.computeDvolume(obj.quadrature);
            xV = obj.quadrature.posgp;
            nGaus  = obj.quadrature.ngaus;
            nF     = size(fVals,1);
            nElem  = size(obj.mesh.connec,1);
            rhs = zeros(nElem,nF);

            % Separate in two loops
            for igaus = 1:nGaus
                xGaus = xV(:,igaus);
                fGaus = fDisc.evaluate(xGaus);
                dVg(:,1) = dV(igaus,:);
                for iF = 1:nF
                    fGausF = squeeze(fGaus(iF,:,:));
                    Ni = 1;
                    int = Ni*fGausF.*dVg;
                    rhs(:,iF) = rhs(:,iF) + int;
                end
            end
        end

    end
end

