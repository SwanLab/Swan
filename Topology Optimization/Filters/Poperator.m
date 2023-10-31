classdef Poperator < handle
    
   properties (Access = public)
       value 
   end
    
   properties (Access = private)
       mesh
       M
   end
    
   methods (Access = public)
       
       function obj = Poperator(cParams)
           obj.init(cParams);
           obj.createMassMatrix();
           obj.createOperator();
       end
       
   end
   
   methods (Access = private)
       
        function init(obj,cParams)
            obj.mesh = cParams.mesh;
        end
       
        function createMassMatrix(obj)
            s.type  = 'MassMatrix';
            s.mesh  = obj.mesh;
            s.test  = P1Function.create(obj.mesh, 1);
            s.trial = P1Function.create(obj.mesh, 1);
            s.quadratureOrder = 'QUADRATICMASS';
            LHS   = LHSintegrator.create(s);
            obj.M = LHS.compute();
        end
       
        function createOperator(obj)
            nelem     = obj.mesh.nelem;
            nodesElem = obj.mesh.nnodeElem;
            ndof      = obj.mesh.nnodes;
            connec    = obj.mesh.connec;
            T = sparse(nelem,ndof);
            for inode = 1:nodesElem
                nodes(:,1) = connec(:,inode);
                I          = ones(nelem,1);
                incT       = sparse(1:nelem,nodes,I,nelem,ndof);
                T          = T + incT;
            end
            m         = T*sum(obj.M,2);
            mInv      = spdiags(1./m,0,length(m),length(m));
            obj.value = mInv*T;
        end

   end

end