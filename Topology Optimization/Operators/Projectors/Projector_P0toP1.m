classdef Projector_P0toP1 < handle
    
    % Eventually should become a single projector

    properties (Access = public)
        value
    end

    properties (Access = private)
        % From Poperator
        mesh
        connec
        nelem
        nnode
        npnod

        field, M
        fieldP0
    end
    
    methods (Access = public)

        function obj = Projector_P0toP1(cParams)
            obj.init(cParams);
            obj.createP1Field();
            obj.createP1MassMatrix();
            obj.createP0Field();
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
    
        %% From Poperator.m
        function createP1Field(obj)
            s.mesh               = obj.mesh;
            s.ndimf              = 3;
            s.interpolationOrder = 'LINEAR';
            s.quadratureOrder    = 'QUADRATICMASS';
            obj.field = Field(s);
        end
       
        function createP1MassMatrix(obj)
            s.type  = 'MassMatrix';
            s.mesh  = obj.mesh;
            s.field = obj.field;
            LHS = LHSintegrator.create(s);
            obj.M = LHS.compute();
        end

        function createP0Field(obj)
            s.mesh               = obj.mesh;
            s.ndimf              = 1;
            s.interpolationOrder = 'CONSTANT';
            f = Field(s);
            obj.fieldP0 = f;
        end

        function rhs = createRHS(obj, fefun)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('LINEAR');
            dV = obj.mesh.computeDvolume(quad);
            xV = quad.posgp;

            obj.mesh.interpolation.computeShapeDeriv(xV);
            shape = obj.field.interpolation.shape;
            nGaus  = quad.ngaus;
            nF     = size(fefun.fElem,1);
            nElem  = size(fefun.fElem,3);
            rhs = zeros(nElem,nF);
            for igaus = 1:nGaus
                fGaus = fefun.interpolateFunction(xV(:,igaus));
                dVg(:,1) = dV(igaus,:);
                for iF = 1:nF
                    fGausF = squeeze(fGaus(iF,:,:));
                    Ni = shape(igaus,:);
                    int = Ni*fGausF.*dVg;
                    rhs(:,iF) = rhs(:,iF) + int;
                end
            end
        end
        
    end

end

