classdef Poperator < handle
    
   properties (Access = public)
       value 
   end
    
   properties (Access = private)
       dim
       mesh
       nelem
       nnode
       npnod
       connec
       diffReacProb
       M
   end
    
   methods (Access = public)
       
       function obj = Poperator(cParams)
           obj.init(cParams);
           obj.computeDimensions();
           obj.createMassMatrix();
           obj.createOperator();
       end
       
   end
   
   methods (Access = private)
       
        function init(obj,cParams)
            obj.nelem  = cParams.nelem;
            obj.nnode  = cParams.nnode;
            obj.npnod  = cParams.npnod;
            obj.connec = cParams.connec;
            obj.mesh   = cParams.diffReactEq.mesh;
        end
       
        function createMassMatrix(obj)
            s.type         = 'MassMatrix';
            s.quadType     = 'QUADRATICMASS';
            s.mesh         = obj.mesh;
            s.globalConnec = obj.mesh.connec;
            s.dim          = obj.dim;
            LHS = LHSintegrator.create(s);
            obj.M = LHS.compute();
        end

        function computeDimensions(obj)
            s.ngaus = [];
            s.mesh  = obj.mesh;
            s.pdim  = '1D';
            d       = DimensionVariables(s);
            d.compute();
            obj.dim = d;
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