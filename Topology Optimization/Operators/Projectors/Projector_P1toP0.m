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
    end
    
    methods (Access = public)

        function obj = Projector_P1toP0(cParams)
            obj.init(cParams);
            obj.createField();
            obj.createOperator();
        end

        function xProj = project(obj, x)
            ndimf  = obj.field.dim.ndimf;
            RHS = zeros(obj.nelem*ndimf,1);
            ngaus = size(x,2);
            for igaus = 1:ngaus
                dvolu = obj.field.geometry.dvolu(:,igaus);
                dv = repmat(dvolu, [ndimf, 1]);
                RHS = RHS + dv.*x(:,igaus); % Why dv here instead of obj.value
            end
            xProj = obj.value*RHS;
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
            s.ndimf              = 3;
            s.interpolationOrder = 'LINEAR';
            s.quadratureOrder    = 'QUADRATICMASS';
            obj.field = Field(s);
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
            obj.value = eye(nelem*ndimf)*T; % Review
        end
        
    end

end

