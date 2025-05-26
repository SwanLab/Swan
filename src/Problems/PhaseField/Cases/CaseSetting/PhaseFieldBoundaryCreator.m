classdef PhaseFieldBoundaryCreator < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        createBoundaryConditions
        boundaryConditions
    end
    
    properties (SetAccess = private, GetAccess = public)
        type
        mesh
        bcValues
        step
    end
    
    methods (Access = public)
        
        function obj = PhaseFieldBoundaryCreator(mesh,cParams)
            obj.init(mesh,cParams)
            obj.defineBoundaryConditions();
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
        
    end
    
    methods (Access = private)
        
        function init(obj,mesh,cParams)
            obj.mesh = mesh;
            obj.type = cParams.type.bc;
            obj.bcValues = cParams.bcValues;
            obj.step = 0;
        end
        
        function defineBoundaryConditions(obj)
           switch obj.type
               case 'bending'
                   obj.createBoundaryConditions = @obj.createBendingConditions;
               case 'forceTraction'
                   obj.createBoundaryConditions = @obj.createForceTractionConditions;
               case 'displacementTraction'
                   obj.createBoundaryConditions = @obj.createDisplacementTractionConditions;
               case 'displacementShear'
                   obj.createBoundaryConditions = @obj.createDisplacementShearConditions;
               case 'displacementMixed'
                   obj.createBoundaryConditions = @obj.createDisplacementMixedConditions;
               case 'Lshape'
                   % obj.createBoundaryConditions = @obj.createLshapeDisplacementConditions;
               case 'FiberMatrix'
                   obj.createBoundaryConditions = @obj.createFiberMatrixDisplacementConditions;
           end
       end


       function createBendingConditions(obj,bcVal)
           isInLeft = @(coord) (abs(coord(:,1) - min(coord(:,1)))< 1e-12);
           sDir.domain    = @(coor) isInLeft(coor);
           sDir.direction = [1,2];
           sDir.value     = 0;
           Dir1 = DirichletCondition(obj.mesh,sDir);

           isInRight = @(coord) (abs(coord(:,1) - max(coord(:,1)))< 1e-12);
           sDir.domain    = @(coor) isInRight(coor);
           sDir.direction = [2];
           sDir.value     = bcVal;
           Dir2 = DirichletCondition(obj.mesh,sDir);

           s.mesh = obj.mesh;
           s.dirichletFun = [Dir1 Dir2];
           s.pointloadFun = [];
           s.periodicFun = [];
           obj.boundaryConditions = BoundaryConditions(s);
        end
        
         function createForceTractionConditions(obj,bcVal)
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

             sNeum.domain    = @(coor) isInUp(coor);
             sNeum.direction = [2];
             sNeum.value     = bcVal;
             Neum1 = PointLoad(obj.mesh,sNeum);

             s.mesh = obj.mesh;
             s.dirichletFun = [Dir1, Dir2];
             s.pointloadFun = [Neum1];
             s.periodicFun = [];
             obj.boundaryConditions = BoundaryConditions(s);           
         end
        
         function createDisplacementTractionConditions(obj,bcVal)
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
             sDir.value     = bcVal;       
             Dir3 = DirichletCondition(obj.mesh,sDir);

             s.mesh = obj.mesh;
             s.dirichletFun = [Dir1 Dir2 Dir3];
             s.pointloadFun = [];
             s.periodicFun = [];
             obj.boundaryConditions = BoundaryConditions(s);
         end

         function createDisplacementShearConditions(obj,bcVal)
             isInDown = @(coor) (abs(coor(:,2) - min(coor(:,2)))  < 1e-12);
             sDir.domain    = @(coor) isInDown(coor);
             sDir.direction = [1,2];
             sDir.value     = 0;
             Dir1 = DirichletCondition(obj.mesh,sDir);

             isInUp = @(coor) (abs(coor(:,2) - max(coor(:,2)))  < 1e-12);
             sDir.domain    = @(coor) isInUp(coor);
             sDir.direction = [1];
             sDir.value     = bcVal;
             Dir2 = DirichletCondition(obj.mesh,sDir);

             sDir.domain    = @(coor) isInUp(coor);
             sDir.direction = [2];
             sDir.value     = 0;
             Dir3 = DirichletCondition(obj.mesh,sDir);             

             s.mesh = obj.mesh;
             s.dirichletFun = [Dir1 Dir2 Dir3];
             s.pointloadFun = [];
             s.periodicFun = [];
             obj.boundaryConditions = BoundaryConditions(s);
         end

         function createDisplacementMixedConditions(obj,bcVal)
             angle = pi/4;

             isInDown = @(coor) (abs(coor(:,2) - min(coor(:,2)))  < 1e-12);
             sDir.domain    = @(coor) isInDown(coor);
             sDir.direction = [1,2];
             sDir.value     = 0;
             Dir1 = DirichletCondition(obj.mesh,sDir);

             isInUp = @(coor) (abs(coor(:,2) - max(coor(:,2)))  < 1e-12);
             sDir.domain    = @(coor) isInUp(coor);
             sDir.direction = [1];
             sDir.value     = bcVal*cos(angle);
             Dir2 = DirichletCondition(obj.mesh,sDir);

             sDir.domain    = @(coor) isInUp(coor);
             sDir.direction = [2];
             sDir.value     = bcVal*sin(angle);
             Dir3 = DirichletCondition(obj.mesh,sDir);             

             s.mesh = obj.mesh;
             s.dirichletFun = [Dir1 Dir2 Dir3];
             s.pointloadFun = [];
             s.periodicFun = [];
             obj.boundaryConditions = BoundaryConditions(s);
         end

         function createLshapeDisplacementConditions(obj,bcVal)
             isInDown = @(coor) (abs(coor(:,2) - min(coor(:,2)))  < 1e-12);
             sDir.domain    = @(coor) isInDown(coor);
             sDir.direction = [1,2];
             sDir.value     = 0;
             Dir1 = DirichletCondition(obj.mesh,sDir);

             isInTip = @(coor) (abs(coor(:,2)-(max(coor(:,2))+min(coor(:,2)))/2) < 1e-12) & (abs(coor(:,1)-max(coor(:,1))) < 30);
             sDir.domain    = @(coor) isInTip(coor);
             sDir.direction = [1];
             sDir.value     = 0;
             Dir2 = DirichletCondition(obj.mesh,sDir);

             sDir.domain    = @(coor) isInTip(coor);
             sDir.direction = [2];
             sDir.value     = bcVal;
             Dir3 = DirichletCondition(obj.mesh,sDir);

             % Merge
             s.mesh = obj.mesh;
             s.dirichletFun = [Dir1 Dir2 Dir3];
             s.pointloadFun = [];
             s.periodicFun = [];
             obj.boundaryConditions = BoundaryConditions(s);
         end

         function createFiberMatrixDisplacementConditions(obj,bcVal)
             % Enforce fixed Dirichlet conditions to the down nodes
             isDown = @(coord) (abs(coord(:,2) - min(coord(:,2)))< 1e-12);
             sDir.domain    = @(coor) isDown(coor);
             sDir.direction = [1,2];
             sDir.value     = 0;
             Dir1 = DirichletCondition(obj.mesh,sDir);

             % Enforce fixed Dirichlet conditions to the fiber nodes
             isFiber = @(coord) ((coord(:,1)-0.5).^2 + (coord(:,2)-0.5).^2) < (0.25^2 + 1e-5);
             sDir.domain    = @(coor) isFiber(coor);
             sDir.direction = [1,2];
             sDir.value     = 0;
             Dir2 = DirichletCondition(obj.mesh,sDir);

             % Enforce roller Dirichlet conditions to the top nodes
             isUp = @(coord) (abs(coord(:,2) - max(coord(:,2)))< 1e-12);
             sDir.domain    = @(coor) isUp(coor);
             sDir.direction = [1];
             sDir.value     = 0;
             Dir3 = DirichletCondition(obj.mesh,sDir);

             % Enforce displacement at the top
             isUp = @(coord) (abs(coord(:,2) - max(coord(:,2)))< 1e-12);
             sDir.domain    = @(coor) isUp(coor);
             sDir.direction = [2];
             sDir.value     = bcVal;
             Dir4 = DirichletCondition(obj.mesh,sDir);

             % Merge
             s.mesh = obj.mesh;
             s.dirichletFun = [Dir1 Dir2 Dir3 Dir4];
             s.pointloadFun = [];
             s.periodicFun = [];
             obj.boundaryConditions = BoundaryConditions(s);
         end

    end
    
end