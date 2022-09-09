classdef Projector_P1toP0 < handle

    properties (Access = public)
        value
    end

    properties (Access = private)
        mesh
        connec
        nelem
        nnode
        npnod
        field
        M
    end

    methods (Access = public)

        function obj = Projector_P1toP0(cParams)
            obj.init(cParams);
            obj.createField();
            obj.createMassMatrix();
        end

        function xProj = project(obj, x)
            RHS = obj.createRHS(x);
            xProj = obj.M\RHS;
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.mesh   = cParams.mesh;
            obj.connec = cParams.connec;
            obj.nelem  = cParams.nelem;
            obj.nnode  = cParams.nnode;
            obj.npnod  = cParams.npnod;
        end

        function createField(obj)
            s.mesh               = obj.mesh;
            s.ndimf              = 1; % ??
            s.interpolationOrder = 'LINEAR';
            s.quadratureOrder    = 'QUADRATICMASS';
            obj.field = Field(s);
        end

        function createMassMatrix(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('LINEAR');
            dv = obj.mesh.computeDvolume(quad);
            obj.M = diag(dv);
        end

        function rhs = createRHS(obj, fefun)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('LINEAR');
            dV = obj.mesh.computeDvolume(quad);
            xV = quad.posgp;

            obj.mesh.interpolation.computeShapeDeriv(xV);
            nGaus  = quad.ngaus;
            nF     = size(fefun.fElem,1);
            nElem  = size(fefun.fElem,3);
            rhs = zeros(nElem,nF);
            for igaus = 1:nGaus
                fGaus = fefun.interpolateFunction(xV(:,igaus));
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

