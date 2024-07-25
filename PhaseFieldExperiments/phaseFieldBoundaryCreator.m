classdef phaseFieldBoundaryCreator < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        createBoundaryConditions
        boundaryConditions
    end
    
    properties (SetAccess = private, GetAccess = public)
        mesh
        bcValues
        step
    end
    
    methods (Access = public)
        
        function obj = phaseFieldBoundaryCreator(mesh,cParams)
            obj.init(mesh,cParams)
            obj.defineBoundaryConditions(cParams);
        end

        function bC = nextStep(obj)
            newStep = obj.step + 1;
            if newStep > length(obj.bcValues)
                disp('Maximum step reached. Previous BC are being used!!!')
            else
                obj.step = newStep;
            end
            prescribedVal = obj.bcValues(obj.step);
            obj.createBoundaryConditions(prescribedVal)
            bC = obj.boundaryConditions;
        end

        function rBC = returnReactionsBoundaryConditions()
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,mesh,cParams)
            obj.mesh = mesh;
            obj.bcValues = cParams.bcValues;
            obj.step = 0;
        end
        
        function defineBoundaryConditions(obj,cParams)
           switch cParams.type.bc
               case 'bending'
                   obj.createBoundaryConditions = @obj.createBendingConditions;
               case 'forceTraction'
                   obj.createBoundaryConditions = @obj.createForceTractionConditions;
               case 'displacementTraction'
                   obj.createBoundaryConditions = @obj.createDisplacementTractionConditions;
               case 'Lshape'
                   obj.createBoundaryConditions = @obj.createLshapeDisplacementConditions;
               case 'FiberMatrix'
                   % obj.createFiberMatrixDisplacementConditions(prescribedVal);
           end
       end


       function createBendingConditions(obj,uVal)
           isInLeft = @(coord) (abs(coord(:,1) - min(coord(:,1)))< 1e-12);
           sDir.domain    = @(coor) isInLeft(coor);
           sDir.direction = [1,2];
           sDir.value     = 0;
           Dir1 = DirichletCondition(obj.mesh,sDir);

           isInRight = @(coord) (abs(coord(:,1) - max(coord(:,1)))< 1e-12);
           sDir.domain    = @(coor) isInRight(coor);
           sDir.direction = [2];
           sDir.value     = uVal;
           Dir2 = DirichletCondition(obj.mesh,sDir);

           s.mesh = obj.mesh;
           s.dirichletFun = [Dir1 Dir2];
           s.pointloadFun = [];
           s.periodicFun = [];
           obj.boundaryConditions = BoundaryConditions(s);
        end
        
         function createForceTractionConditions(obj,fVal)
             isInDown = @(coor) (abs(coor(:,2) - min(coor(:,2)))  < 1e-12);
             sDir.domain    = @(coor) isInDown(coor);
             sDir.direction = [1,2];
             sDir.value     = 0;
             Dir1 = DirichletCondition(obj.mesh,sDir);

             isInUp = @(coor) (abs(coor(:,2) - max(coor(:,2)))  < 1e-12);
             sNeum.domain    = @(coor) isInUp(coor);
             sNeum.direction = [2];
             sNeum.value     = fVal;
             Neum1 = PointLoad(obj.mesh,sNeum);

             s.mesh = obj.mesh;
             s.dirichletFun = [Dir1];
             s.pointloadFun = [Neum1];
             s.periodicFun = [];
             obj.boundaryConditions = BoundaryConditions(s);           
         end
        
         function createDisplacementTractionConditions(obj,uVal)
             isInDown = @(coor) (abs(coor(:,2) - min(coor(:,2)))  < 1e-12);
             sDir.domain    = @(coor) isInDown(coor);
             sDir.direction = [1,2];
             sDir.value     = 0;
             Dir1 = DirichletCondition(obj.mesh,sDir);
             
             isInUp = @(coor) (abs(coor(:,2) - max(coor(:,2)))  < 1e-12);
             sDir.domain    = @(coor) isInUp(coor);
             sDir.direction = [1];
             sDir.value     = 0;
             Dir2 = DirichletCondition(obj.mesh,sDir);

             sDir.domain    = @(coor) isInUp(coor);
             sDir.direction = [2];
             sDir.value     = uVal;       
             Dir3 = DirichletCondition(obj.mesh,sDir);

             s.mesh = obj.mesh;
             s.dirichletFun = [Dir1 Dir2 Dir3];
             s.pointloadFun = [];
             s.periodicFun = [];
             obj.boundaryConditions = BoundaryConditions(s);
         end
         
         function createLshapeDisplacementConditions(obj,uVal)
             isInDown = @(coor) (abs(coor(:,2) - min(coor(:,2)))  < 1e-12);
             sDir.domain    = @(coor) isInDown(coor);
             sDir.direction = [1,2];
             sDir.value     = 0;
             Dir1 = DirichletCondition(obj.mesh,sDir);

             % Enforce displacement under tip
             isInTip = @(coor) (abs(coor(:,2)-(max(coor(:,2))+min(coor(:,2)))/2) < 1e-12) & (abs(coor(:,1)-max(coor(:,1))) < 1e-12);

             sDir.domain    = @(coor) isInTip(coor);
             sDir.direction = [2];
             sDir.value     = uVal;
             Dir2 = DirichletCondition(obj.mesh,sDir);

             % Merge
             s.mesh = obj.mesh;
             s.dirichletFun = [Dir1 Dir2];
             s.pointloadFun = [];
             s.periodicFun = [];
             obj.boundaryConditions = BoundaryConditions(s);
         end

         function createFiberMatrixDisplacementConditions(obj,uVal)
             % nodes = 1:obj.mesh.nnodes;
             % ndim = 2;
             % 
             % % Enforce fixed Dirichlet conditions to the down nodes
             % downSide = min(obj.mesh.coord(:,2));
             % isInDown = abs(obj.mesh.coord(:,2)-downSide) < 1e-12;
             % dirichletDown = zeros(ndim*length(nodes(isInDown)),3);
             % for i=1:ndim
             %     dirichletDown(i:2:end,1) = nodes(isInDown);
             %     dirichletDown(i:2:end,2) = i;
             % end
             % 
             % % Enforce fixed Dirichlet conditions to the fiber nodes
             % center = [(min(obj.mesh.coord(:,1))+max(obj.mesh.coord(:,1)))/2;
             %           (min(obj.mesh.coord(:,2))+max(obj.mesh.coord(:,2)))/2];
             % radius = 0.2;
             % isInCircle = ((obj.mesh.coord(:,1)-center(1)).^2 + (obj.mesh.coord(:,2)-center(2)).^2) ...
             %              < (radius^2 + 1e-5);
             % dirichletCircle = zeros(ndim*length(nodes(isInCircle)),3);
             % for i=1:ndim
             %     dirichletCircle(i:2:end,1) = nodes(isInCircle);
             %     dirichletCircle(i:2:end,2) = i;
             % end
             % 
             % % Enforce roller Dirichlet conditions to the top nodes
             % upSide  = max(obj.mesh.coord(:,2));
             % isInUp = abs(obj.mesh.coord(:,2)-upSide)< 1e-12;
             % dirichletUp   = zeros(ndim*length(nodes(isInUp)),3);
             % for i=1:ndim
             %     dirichletUp(i:2:end,1) = nodes(isInUp);
             %     dirichletUp(i:2:end,2) = i;
             % end
             % 
             % % Enforce displacement at the top
             % dirichletUp(2:2:end,3) = uVal;
             % 
             % % Merge
             % bc.dirichlet = [dirichletDown; dirichletUp; dirichletCircle];
             % bc.pointload = [];
             % obj.boundaryConditions = bc;
         end

    end
    
end