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
%             obj.createP0MassMatrix();
            obj.createOperator();
        end

        function xProj = project(obj, x)
            ndimf  = obj.field.dim.ndimf;
            RHS = zeros(obj.nelem*ndimf,1);
            ngaus = size(x,2);
            for igaus = 1:ngaus
                dvolu = obj.field.geometry.dvolu(:,igaus);
                dv = repmat(dvolu, [ndimf, 1]);
                RHS = RHS + dv.*x(:,igaus);
            end
            xProj = obj.value'*RHS;
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

        function computeLHS(obj)
        end

        function computeRHS(obj)
            % ez
        end

        function solve(obj)
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
       
        function createP0MassMatrix(obj)
            s.type  = 'MassMatrix';
            s.mesh  = obj.mesh;
            s.field = obj.fieldP0;
            LHS = LHSintegrator.create(s);
            obj.M = LHS.compute();
        end
       
        function createOperator(obj)
            nelem  = obj.nelem;
            nnode  = obj.nnode;
            npnod  = obj.npnod;
            ndimf  = obj.field.dim.ndimf;
            connec = obj.connec;
            T = sparse(nelem*ndimf,npnod*ndimf);
            for inode = 1:nnode
                dofsArr = [];
                nods = connec(:, inode);
                for idim = 1:ndimf
                    dofs  = ndimf*(nods - 1) + idim;
                    dofsArr = [dofsArr; dofs]; 
                end
                nodes(:,1) = dofsArr;
                I = ones(nelem*ndimf,1);
                incT = sparse(1:nelem*ndimf,nodes,I,nelem*ndimf,npnod*ndimf);
                T = T + incT;
            end
            m = T*sum(obj.M,2);
            mInv = spdiags(1./m,0,length(m),length(m));
            obj.value = mInv*T;
        end
        
    end

end

