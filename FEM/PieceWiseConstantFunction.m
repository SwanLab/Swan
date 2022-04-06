classdef PieceWiseConstantFunction < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        integrator
        quadrature
        dim
        dofsInElem
    end
    
    properties (Access = private)
       mesh 
       fValues
       quadOrder
    end
    
    methods (Access = public)
        
        function obj = PieceWiseConstantFunction(cParams)
            obj.init(cParams);
            obj.createQuadrature();
            obj.createDimensions();
            obj.computeDofConnectivity();
        end
        
        function fNodal = projectToLinearNodalFunction(obj)
            obj.createIntegrator();
            LHS = obj.computeLHS();
            RHS = obj.computeRHS();
            fNodal = (LHS\RHS);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh      = cParams.mesh;
            obj.fValues   = cParams.fValues;
            obj.quadOrder = 'LINEAR';
        end
        
        function q = createQuadrature(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature(obj.quadOrder);
            obj.quadrature = q;
        end

        function createDimensions(obj)
            q = Quadrature();
            q = q.set(obj.mesh.type);
            s.mesh = obj.mesh;
            s.pdim = '1D';
            s.ngaus = q.ngaus;
            d = DimensionVariables(s);
            d.compute();
            obj.dim = d;
        end

        function createIntegrator(obj)
            s.type = 'SIMPLE';
            s.mesh = obj.mesh;
            s.npnod = obj.mesh.npnod;
            s.globalConnec = obj.mesh.connec;
            int = Integrator.create(s);
            obj.integrator = int;
        end
        
        function x = computeXgauss(obj)
            xG = obj.quadrature.posgp;
            x = repmat(xG,[1,1,obj.mesh.nelem]);
        end
        
        function f = computeFgauss(obj)
            ngaus = obj.quadrature.ngaus;
            fV(1,:) = obj.fValues;
            f = repmat(fV,[ngaus,1]);
        end
        
        function LHS = computeLHS(obj)
            s.mesh         = obj.mesh;
            s.globalConnec = obj.mesh.connec;
            s.npnod        = obj.mesh.npnod;
            s.type         = 'MassMatrix';
            s.dim          = obj.dim;
            s.quadType     = 'QUADRATIC';
            s.dofsInElem   = obj.dofsInElem;
            lhs = LHSintegrator.create(s);
            LHS = lhs.compute();
        end

       function computeDofConnectivity(obj)
           connec = obj.mesh.connec;
           ndimf  = obj.dim.ndimField;
           nnode  = obj.dim.nnode;
           dofsElem  = zeros(nnode*ndimf,size(connec,1));
           for inode = 1:nnode
               for iunkn = 1:ndimf
                   idofElem   = obj.nod2dof(inode,iunkn);
                   globalNode = connec(:,inode);
                   idofGlobal = obj.nod2dof(globalNode,iunkn);
                   dofsElem(idofElem,:) = idofGlobal;
               end
           end
           obj.dofsInElem = dofsElem;
       end        
        
        function RHS = computeRHS(obj)
            fG = obj.computeFgauss;
            xG = obj.computeXgauss;
            RHS = obj.integrator.integrateFgauss(fG,xG,obj.quadOrder);
        end

        function idof = nod2dof(obj, inode, iunkn)
            ndimf = obj.dim.ndimField;
            idof(:,1)= ndimf*(inode - 1) + iunkn;
        end         
        
    end
    
end