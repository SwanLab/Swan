classdef BoundaryContionsForPhaseFieldCreator < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        boundaryConditions
    end
    
    properties (Access = private)
        mesh
    end
    
    methods (Access = public)
        
        function obj = BoundaryContionsForPhaseFieldCreator(cParams)
            obj.init(cParams)
            
        end
        
        function bC = create(obj,iter,nIter)
            obj.createBoundaryConditions(iter,nIter)
            bC = obj.boundaryConditions;            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh = cParams.mesh;
        end
        
       function createBoundaryConditions(obj,iter,nIter)
          %  obj.createBendingConditions(iter,nIter)
          %  obj.createForceTractionConditions(iter,nIter)
            obj.createDisplacementUniaxialConditions(iter)
        end
        
        
        function createBendingConditions(obj,iter,nIter)
            Force = 1.4;
            dirichletNodes = abs(obj.mesh.coord(:,1)-0) < 1e-12;
            rightSide  = max(obj.mesh.coord(:,1));
            isInRight = abs(obj.mesh.coord(:,1)-rightSide)< 1e-12;
            isInMiddleEdge = abs(obj.mesh.coord(:,2)-1.5) < 0.1;
            forceNodes = isInRight & isInMiddleEdge;
            nodes = 1:obj.mesh.nnodes;
            
            ndim = 2;
            bc.dirichlet = zeros(ndim*length(nodes(dirichletNodes)),3);
            for i=1:ndim
                bc.dirichlet(i:2:end,1) = nodes(dirichletNodes);
                bc.dirichlet(i:2:end,2) = i;
            end
            bc.dirichlet(:,3) = 0;
            
            bc.pointload(:,1) = nodes(forceNodes);
            bc.pointload(:,2) = 2;
            bc.pointload(:,3) = -Force*(iter/nIter);
            obj.boundaryConditions = bc;
        end
        
         function createForceTractionConditions(obj,iter,nIter)
            Force = 1.4;
            DownSide = min(obj.mesh.coord(:,2));
            isInDown = abs(obj.mesh.coord(:,2)-DownSide) < 1e-12;
            UpSide  = max(obj.mesh.coord(:,2));
            isInUp = abs(obj.mesh.coord(:,2)-UpSide)< 1e-12;
            forceNodes = isInUp;
            nodes = 1:obj.mesh.nnodes;
            
            ndim = 2;
            bc.dirichlet = zeros(ndim*length(nodes(isInDown)),3);
            for i=1:ndim
                bc.dirichlet(i:2:end,1) = nodes(isInDown);
                bc.dirichlet(i:2:end,2) = 1;
            end
            
            bc.dirichlet(2:2:end,3) = 0;
            bc.dirichlet = [1 1 0; 
                            1 2 0;
                            2 1 0;
                            2 2 0;
                            3 1 0;
                            4 1 0];                     
            
            bc.pointload(:,1) = nodes(forceNodes);
            bc.pointload(:,2) = 2;
            bc.pointload(:,3) = Force*(iter/nIter);
            obj.boundaryConditions = bc;            
         end
        
         function createDisplacementUniaxialConditions(obj,iter)
            DisplacementStep = 1e-4;
            DownSide = min(obj.mesh.coord(:,2));
            isInDown = abs(obj.mesh.coord(:,2)-DownSide) < 1e-12;
            UpSide  = max(obj.mesh.coord(:,2));
            isInUp = abs(obj.mesh.coord(:,2)-UpSide)< 1e-12;
            nodes = 1:obj.mesh.nnodes;

            ndim = 2;
            DirichletDown = zeros(ndim*length(nodes(isInDown)),3);
            DirichletUp   = zeros(ndim*length(nodes(isInUp)),3);
            for i=1:ndim
                DirichletDown(i:2:end,1) = nodes(isInDown);
                DirichletDown(i:2:end,2) = i;
                DirichletDown(:,3) = 0;

               DirichletUp(i:2:end,1) = nodes(isInUp);
               DirichletUp(i:2:end,2) = i;                
            end            
            DirichletUp(2:2:end,3) = DisplacementStep*(iter);

            bc.dirichlet = [DirichletDown; DirichletUp];
            bc.pointload = [];

            obj.boundaryConditions = bc;             
         end        
        
    end
    
end