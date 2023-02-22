classdef MinimumDiscGradFieldWithVectorInH1 < handle

    properties (Access = private)
        LHS
        RHS
    end
    
    properties (Access = private)
        mesh
        rhsType
        meshCont
        fValue
        interp
        field
        interpolator
    end
    
    methods (Access = public)
        
        function obj = MinimumDiscGradFieldWithVectorInH1(cParams)
            obj.init(cParams);
            obj.createField();
        end

        function u = solve(obj)
            obj.computeLHS();
            obj.computeRHS();
            uC = obj.solveSystem();
            In = obj.interpolator; 
            u  = In*uC; 
            u = reshape(u,obj.mesh.nnodeElem,[])'; % Eh          
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.meshCont     = cParams.mesh;
            obj.rhsType      = cParams.rhsType;
            obj.mesh         = obj.meshCont.createDiscontinuousMesh();            
            obj.fValue       = cParams.fValue;
            obj.interpolator = cParams.interpolator;
        end

        function createField(obj)
            s.mesh               = obj.mesh;
            s.ndimf              = 1;
            s.interpolationOrder = obj.mesh.interpolation.order;
            obj.field = Field(s);
        end
%         
        function computeLHS(obj)
            K = obj.computeStiffnessMatrix();
            In = obj.interpolator;
            K = In'*K*In;
            obj.LHS = K;
        end
        
        function K = computeStiffnessMatrix(obj)
            s.mesh         = obj.mesh;
            s.globalConnec = obj.mesh.connec;
            s.type         = 'StiffnessMatrix';
            s.field        = obj.field;

            lhs = LHSintegrator.create(s);
            K = lhs.compute();
        end

        function fGauss = computeFGauss(obj)%%%Ehhhhh
            cV = obj.fValue;
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('QUADRATIC');
            int = Interpolation.create(obj.mesh,obj.mesh.interpolation.order);
            int.computeShapeDeriv(q.posgp);
            s.mesh = obj.mesh;
            g = Geometry.create(s);
            g.computeGeometry(q,int);
            dN = g.dNdx;

            %dN = int.deriv;
            nDim = size(dN,1);
            nnode = size(dN,2);
            fGauss = zeros(nDim,q.ngaus,obj.mesh.nelem);
            for igaus = 1:q.ngaus
                for idim = 1:nDim
                    for iNode = 1:nnode
                        dNi = squeeze(dN(idim,iNode,:,igaus));
                        grad = dNi.*squeeze(cV(1,iNode,:));
                        fG(:,1) = squeeze(fGauss(idim,igaus,:));
                        fGauss(idim,igaus,:) = fG + grad;
                    end
                end
            end

        end
        
        function computeRHS(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('QUADRATIC');
            s.fType     = 'Gauss';
            s.fGauss    = obj.computeFGauss();
            s.xGauss    = q.posgp;
            s.mesh      = obj.mesh;
            s.type      = obj.mesh.type;
            s.quadOrder = q.order;
            s.npnod     = obj.field.dim.ndofs;
            s.type      = obj.rhsType;
            s.globalConnec = obj.mesh.connec;
            rhs  = RHSintegrator.create(s);
            rhsV = rhs.compute();        
            In = obj.interpolator;
            rhsV = In'*rhsV;
            %obj.RHS = [rhsV;0];
            obj.RHS = rhsV;
        end
        
        function f = assembleIntegrand(obj,rhsCells)
            integrand = rhsCells;
            ndofs  = obj.mesh.nnodes;
            connec = obj.mesh.connec;
            nnode  = size(connec,2);
            f = zeros(ndofs,1);
            for inode = 1:nnode
                int = integrand(:,inode);
                con = connec(:,inode);
                f = f + accumarray(con,int,[ndofs,1],@sum,0);
            end
        end
        
        function u = solveSystem(obj)
            a.type = 'DIRECT';
            s = Solver.create(a);
            u = s.solve(obj.LHS,obj.RHS);
            u = u(1:end);
        end
        
    end
    
end