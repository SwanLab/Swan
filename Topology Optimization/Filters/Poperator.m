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
       field
   end
    
   methods (Access = public)
       
       function obj = Poperator(cParams)
           obj.init(cParams);
           obj.createField();
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