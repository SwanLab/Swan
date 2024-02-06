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
        
        function bC = create(obj,prescribedVal)
            obj.createBoundaryConditions(prescribedVal)
            bC = obj.boundaryConditions;            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh = cParams.mesh;
        end
        
       function createBoundaryConditions(obj,prescribedVal)
           obj.createBendingConditions(prescribedVal)
          % obj.createForceTractionConditions(prescribedVal)
          % obj.createDisplacementTractionConditions(prescribedVal)
          % obj.createLshapeDisplacementConditions(prescribedVal)
          % obj.createFiberMatrixDisplacementConditions(prescribedVal);
        end
        

        function createBendingConditions(obj,uVal)
            nodes = 1:obj.mesh.nnodes;
            ndim = 2;

            % Enforce fixed Dirichlet conditions to the right nodes
            leftSide  = min(obj.mesh.coord(:,1));
            isInLeft = abs(obj.mesh.coord(:,1)-leftSide)< 1e-12;
            dirichletLeft     = zeros(ndim*length(nodes(isInLeft)),3); 
            for i=1:ndim
                dirichletLeft(i:2:end,1) = nodes(isInLeft);
                dirichletLeft(i:2:end,2) = i;
            end

            % Enforce displacement at the tip
            rightSide = max(obj.mesh.coord(:,1));
            isInRight = abs(obj.mesh.coord(:,1)-rightSide) < 1e-12;
            dirichletRight = zeros(length(nodes(isInRight)),3);
            dirichletRight(:,1) = nodes(isInRight);
            dirichletRight(:,2) = 2;
            dirichletRight(:,3) = uVal;
            
            % Merge all B.C.
            bc.dirichlet = [dirichletLeft; dirichletRight];
            bc.pointload = [];
            obj.boundaryConditions = bc;
        end
        
         function createForceTractionConditions(obj,fVal)
            nodes = 1:obj.mesh.nnodes;
            ndim = 2; 

            % Enforce fixed Dirichlet conditions to the down nodes
            downSide = min(obj.mesh.coord(:,2));
            isInDown = abs(obj.mesh.coord(:,2)-downSide) < 1e-12;
            dirichletDown = zeros(ndim*length(nodes(isInDown)),3);
            for i=1:ndim
                dirichletDown(i:2:end,1) = nodes(isInDown);
                dirichletDown(i:2:end,2) = i;
            end

            % Enforce roller Dirichlet conditions to the top nodes
            upSide  = max(obj.mesh.coord(:,2));
            isInUp = abs(obj.mesh.coord(:,2)-upSide)< 1e-12;
            dirichletUp   = zeros(length(nodes(isInUp)),3);
            dirichletUp(:,1) = nodes(isInUp);
            dirichletUp(:,2) = 1;

            % Enforce force at the top
            bc.pointload(:,1) = nodes(isInUp);
            bc.pointload(:,2) = 2;
            bc.pointload(:,3) = fVal/length(nodes(isInUp));

            % Merge
            bc.dirichlet = [dirichletDown; dirichletUp];                  
            obj.boundaryConditions = bc;            
         end
        
         function createDisplacementTractionConditions(obj,uVal)
             nodes = 1:obj.mesh.nnodes;
             ndim = 2;

             % Enforce fixed Dirichlet conditions to the down nodes
             downSide = min(obj.mesh.coord(:,2));
             isInDown = abs(obj.mesh.coord(:,2)-downSide) < 1e-12;
             dirichletDown = zeros(ndim*length(nodes(isInDown)),3);
             for i=1:ndim
                 dirichletDown(i:2:end,1) = nodes(isInDown);
                 dirichletDown(i:2:end,2) = i;
             end

             % Enforce roller Dirichlet conditions to the top nodes
             upSide  = max(obj.mesh.coord(:,2));
             isInUp = abs(obj.mesh.coord(:,2)-upSide)< 1e-12;
             dirichletUp   = zeros(ndim*length(nodes(isInUp)),3);
             for i=1:ndim
                 dirichletUp(i:2:end,1) = nodes(isInUp);
                 dirichletUp(i:2:end,2) = i;
             end

             % Enforce displacement at the top
             dirichletUp(2:2:end,3) = uVal;

             % Merge
             bc.dirichlet = [dirichletDown; dirichletUp];
             bc.pointload = [];
             obj.boundaryConditions = bc;
         end
         
         function createLshapeDisplacementConditions(obj,uVal)
             nodes = 1:obj.mesh.nnodes;
             ndim = 2;

             % Enforce fixed Dirichlet conditions to the down nodes
             downSide = min(obj.mesh.coord(:,2));
             isInDown = abs(obj.mesh.coord(:,2)-downSide) < 1e-12;
             dirichletDown = zeros(ndim*length(nodes(isInDown)),3);
             for i=1:ndim
                 dirichletDown(i:2:end,1) = nodes(isInDown);
                 dirichletDown(i:2:end,2) = i;
             end

             % Enforce displacement under tip
             rightSide = max(obj.mesh.coord(:,1));
             isInRight = abs(obj.mesh.coord(:,1)-rightSide) < 1e-12;
             upSide = max(obj.mesh.coord(:,2));
             midPlane = upSide - downSide;
             isInMiddle = abs(obj.mesh.coord(:,2)-midPlane) < 1e-12;
             isInTip = isInRight & isInMiddle;

             dirichletTip   = zeros(length(nodes(isInTip)),3);
             dirichletTip(:,1) = nodes(isInTip);
             dirichletTip(:,2) = 2;
             dirichletTip(:,3) = uVal;

             % Merge
             bc.dirichlet = [dirichletDown; dirichletTip];
             bc.pointload = [];
             obj.boundaryConditions = bc;
         end

         function createFiberMatrixDisplacementConditions(obj,uVal)
             nodes = 1:obj.mesh.nnodes;
             ndim = 2;

             % Enforce fixed Dirichlet conditions to the down nodes
             downSide = min(obj.mesh.coord(:,2));
             isInDown = abs(obj.mesh.coord(:,2)-downSide) < 1e-12;
             dirichletDown = zeros(ndim*length(nodes(isInDown)),3);
             for i=1:ndim
                 dirichletDown(i:2:end,1) = nodes(isInDown);
                 dirichletDown(i:2:end,2) = i;
             end

             % Enforce fixed Dirichlet conditions to the fiber nodes
             center = [(min(obj.mesh.coord(:,1))+max(obj.mesh.coord(:,1)))/2;
                       (min(obj.mesh.coord(:,2))+max(obj.mesh.coord(:,2)))/2];
             radius = 0.2;
             isInCircle = ((obj.mesh.coord(:,1)-center(1)).^2 + (obj.mesh.coord(:,2)-center(2)).^2) ...
                          < (radius^2 + 1e-5);
             dirichletCircle = zeros(ndim*length(nodes(isInCircle)),3);
             for i=1:ndim
                 dirichletCircle(i:2:end,1) = nodes(isInCircle);
                 dirichletCircle(i:2:end,2) = i;
             end

             % Enforce roller Dirichlet conditions to the top nodes
             upSide  = max(obj.mesh.coord(:,2));
             isInUp = abs(obj.mesh.coord(:,2)-upSide)< 1e-12;
             dirichletUp   = zeros(ndim*length(nodes(isInUp)),3);
             for i=1:ndim
                 dirichletUp(i:2:end,1) = nodes(isInUp);
                 dirichletUp(i:2:end,2) = i;
             end

             % Enforce displacement at the top
             dirichletUp(2:2:end,3) = uVal;

             % Merge
             bc.dirichlet = [dirichletDown; dirichletUp; dirichletCircle];
             bc.pointload = [];
             obj.boundaryConditions = bc;
         end

    end
    
end