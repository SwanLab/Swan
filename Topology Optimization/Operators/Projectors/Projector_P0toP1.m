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

        field, M,
    end
    
    methods (Access = public)

        function obj = Projector_P0toP1(cParams)
            obj.init(cParams);
            obj.createField();
            obj.createMassMatrix();
            obj.createOperator();
        end

        function xProj = project(obj, x)
            RHS = zeros(obj.nelem,1);
            ng = size(x,2);
            for igaus = 1:ng
                dvolu = obj.field.geometry.dvolu(:,igaus);
                RHS = RHS + dvolu.*x(:,igaus);
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
        function createField(obj)
            s.mesh               = obj.mesh;
            s.ndimf              = 1;
            s.interpolationOrder = 'LINEAR';
            s.quadratureOrder    = 'QUADRATICMASS';
            obj.field = Field(s);
        end
       
        function createMassMatrix(obj)
            s.type  = 'MassMatrix';
            s.mesh  = obj.mesh;
            s.field = obj.field;
            LHS = LHSintegrator.create(s);
            obj.M = LHS.compute();
        end
       
        function createOperator(obj)
            nelem  = obj.nelem;
            nnode  = obj.nnode;
            npnod  = obj.npnod;
            connec = obj.connec;
            T = sparse(nelem,npnod);
            for inode = 1:nnode
                nodes(:,1) = connec(:,inode);
                I = ones(nelem,1);
                incT = sparse(1:nelem,nodes,I,nelem,npnod);
                T = T + incT;
            end
            m = T*sum(obj.M,2);
            mInv = spdiags(1./m,0,length(m),length(m));
            obj.value = mInv*T;
        end
        
    end

end

