classdef BoundaryConditionsCreator < handle

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

        function obj = BoundaryConditionsCreator(mesh,cParams)
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
            obj.mesh     = mesh;
            obj.type     = cParams.type;
            obj.bcValues = cParams.values;
            obj.step     = 0;
        end

        function defineBoundaryConditions(obj)
            switch obj.type
                case 'ForceTractionX'
                    obj.createBoundaryConditions = @obj.createForceTractionXConditions;
                case 'DisplacementTractionX'
                    obj.createBoundaryConditions = @obj.createDisplacementTractionXConditions;
                case 'ForceTractionY'
                    obj.createBoundaryConditions = @obj.createForceTractionYConditions;
                case 'DisplacementTractionY'
                    obj.createBoundaryConditions = @obj.createDisplacementTractionYConditions;
                case 'ThreePointBending'
                    obj.createBoundaryConditions = @obj.createThreePointBendingConditions;
                case 'ForceCompressionZ'
                    obj.createBoundaryConditions = @obj.createForceTractionZConditions;
                case 'DisplacementShear'
                    obj.createBoundaryConditions = @obj.createDisplacementShearConditions;
                case 'DisplacementMixed'
                    obj.createBoundaryConditions = @obj.createDisplacementMixedConditions;
                case 'ForceBending'
                    obj.createBoundaryConditions = @obj.createForceBendingConditions;
                case 'DisplacementBending'
                    obj.createBoundaryConditions = @obj.createDisplacementBendingConditions;
            end
        end

        function createForceTractionXConditions(obj,fVal)
            isLeft = @(coor)  abs(coor(:,1)-min(coor(:,1))) < 1e-12;
            sDir.domain    = @(coor) isLeft(coor);
            sDir.direction = [1,2];
            sDir.value     = 0;
            Dir1 = DirichletCondition(obj.mesh,sDir);

            isRight = @(coor)  abs(coor(:,1)-max(coor(:,1))) < 1e-12;
            sNeum.domain    = @(coor) isRight(coor);
            sNeum.direction = [1];
            sNeum.value     = fVal;
            Neum1 = PointLoad(obj.mesh,sNeum);
            % Remember change bMesh{2} in extWorkFunctional

            s.mesh         = obj.mesh;
            s.dirichletFun = [Dir1];
            s.pointloadFun = [Neum1];
            s.periodicFun  = [];
            obj.boundaryConditions = BoundaryConditions(s);
        end

        function createDisplacementTractionXConditions(obj,uVal)
            isLeft = @(coor)  abs(coor(:,1)-min(coor(:,1))) < 1e-12;
            sDir.domain    = @(coor) isLeft(coor);
            sDir.direction = [1,2];
            sDir.value     = 0;
            Dir1 = DirichletCondition(obj.mesh,sDir);

            isRight = @(coor)  abs(coor(:,1)-max(coor(:,1))) < 1e-12;
            sDir.domain    = @(coor) isRight(coor);
            sDir.direction = [2];
            sDir.value     = 0;
            Dir2 = DirichletCondition(obj.mesh,sDir);

            sDir.domain    = @(coor) isRight(coor);
            sDir.direction = [1];
            sDir.value     = uVal;
            Dir3 = DirichletCondition(obj.mesh,sDir);

            s.mesh         = obj.mesh;
            s.dirichletFun = [Dir1 Dir2 Dir3];
            s.pointloadFun = [];
            s.periodicFun  = [];
            obj.boundaryConditions = BoundaryConditions(s);
        end

        function createForceTractionYConditions(obj,fVal)
            isDown = @(coor) abs(coor(:,2) - min(coor(:,2))) < 1e-12;
            sDir.domain    = @(coor) isDown(coor);
            sDir.direction = [1,2];
            sDir.value     = 0;
            Dir1 = DirichletCondition(obj.mesh,sDir);

            isUp = @(coor) abs(coor(:,2) - max(coor(:,2))) < 1e-12;
            sNeum.domain    = @(coor) isUp(coor);
            sNeum.direction = [2];
            sNeum.value     = fVal;
            Neum1 = TractionLoad(obj.mesh,sNeum,'DIRAC');
            % Remember change bMesh{4} in extWorkFunctional

            s.mesh = obj.mesh;
            s.dirichletFun = [Dir1];
            s.pointloadFun = [Neum1];
            s.periodicFun  = [];
            obj.boundaryConditions = BoundaryConditions(s);
        end

        function createDisplacementTractionYConditions(obj,uVal)
            isDown = @(coor) abs(coor(:,2) - min(coor(:,2)))  < 1e-12;
            sDir.domain    = @(coor) isDown(coor);
            sDir.direction = [1,2];
            sDir.value     = 0;
            Dir1 = DirichletCondition(obj.mesh,sDir);

            isUp = @(coor) abs(coor(:,2) - max(coor(:,2)))  < 1e-12;
            % sDir.domain    = @(coor) isUp(coor);
            % sDir.direction = [1];
            % sDir.value     = 0;
            % Dir2 = DirichletCondition(obj.mesh,sDir);

            sDir.domain    = @(coor) isUp(coor);
            sDir.direction = [2];
            sDir.value     = uVal;
            Dir3 = DirichletCondition(obj.mesh,sDir);

            s.mesh = obj.mesh;
            s.dirichletFun = [Dir1 Dir3];
            s.pointloadFun = [];
            s.periodicFun = [];
            obj.boundaryConditions = BoundaryConditions(s);
        end

        function createThreePointBendingConditions(obj,fVal)
            isLeft = @(coor)  abs(coor(:,1)-min(coor(:,1))) < 1e-12;
            sDir.domain    = @(coor) isLeft(coor);
            sDir.direction = [1,2];
            sDir.value     = 0;
            Dir1 = DirichletCondition(obj.mesh,sDir);

            isRight = @(coor)  abs(coor(:,1)-max(coor(:,1))) < 1e-12;
            sDir.domain    = @(coor) isRight(coor);
            sDir.direction = [1,2];
            sDir.value     = 0;
            Dir2 = DirichletCondition(obj.mesh,sDir);

            isTop = @(coor) abs(coor(:,2)-max(coor(:,2))) < 1e-12;
            isMiddle = @(coor) abs(coor(:,1)-(min(coor(:,1)) + max(coor(:,1)))/2) < 1e-12;
            sNeum.domain    = @(coor) isTop(coor) & isMiddle(coor);
            sNeum.direction = [2];
            sNeum.value     = fVal;
            Neum1 = DirichletCondition(obj.mesh,sNeum);
            % Remember change bMesh{4} in extWorkFunctional

            s.mesh         = obj.mesh;
            s.dirichletFun = [Dir1 Dir2];
            s.pointloadFun = [Neum1];
            s.periodicFun  = [];
            obj.boundaryConditions = BoundaryConditions(s);
        end

        function createForceTractionZConditions(obj,fVal)
            isBottom = @(coor) coor(:,3)-min(coor(:,3)) < 1e-12;
            sDir.domain    = @(coor) isBottom(coor);
            sDir.direction = [1,2,3];
            sDir.value     = 0;
            Dir1 = DirichletCondition(obj.mesh,sDir);

            isTop = @(coor) coor(:,3)-max(coor(:,3)) < 1e-12;
            isSquare = @(coor) (coor(:,1)-(0.75*max(coor(:,1)+0.25*min(coor(:,1)))) < 1e-12) & ...
                               (coor(:,1)-(0.25*max(coor(:,1)+0.75*min(coor(:,1)))) < 1e-12) & ...
                               (coor(:,2)-(0.75*max(coor(:,2)+0.25*min(coor(:,2)))) < 1e-12) & ...
                               (coor(:,2)-(0.25*max(coor(:,2)+0.75*min(coor(:,2)))) < 1e-12);
            sNeum.domain    = @(coor) isTop(coor) & isSquare(coor);
            sNeum.direction = [3];
            sNeum.value     = fVal;
            Neum1 = DirichletCondition(obj.mesh,sNeum);
            % Remember set bMesh{6} in extWorkFunctional

            s.mesh         = obj.mesh;
            s.dirichletFun = [Dir1];
            s.pointloadFun = [Neum1];
            s.periodicFun  = [];
            obj.boundaryConditions = BoundaryConditions(s);
        end

        function createDisplacementShearConditions(obj,uVal)
             isInDown = @(coor) (abs(coor(:,2) - min(coor(:,2)))  < 1e-12);
             sDir.domain    = @(coor) isInDown(coor);
             sDir.direction = [1,2];
             sDir.value     = 0;
             Dir1 = DirichletCondition(obj.mesh,sDir);

             isInUp = @(coor) (abs(coor(:,2) - max(coor(:,2)))  < 1e-12);
             sDir.domain    = @(coor) isInUp(coor);
             sDir.direction = [1];
             sDir.value     = uVal;
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

         function createDisplacementMixedConditions(obj,uVal)
             angle = pi/4;

             isInDown = @(coor) (abs(coor(:,2) - min(coor(:,2)))  < 1e-12);
             sDir.domain    = @(coor) isInDown(coor);
             sDir.direction = [1,2];
             sDir.value     = 0;
             Dir1 = DirichletCondition(obj.mesh,sDir);

             isInUp = @(coor) (abs(coor(:,2) - max(coor(:,2)))  < 1e-12);
             sDir.domain    = @(coor) isInUp(coor);
             sDir.direction = [1];
             sDir.value     = uVal*cos(angle);
             Dir2 = DirichletCondition(obj.mesh,sDir);

             sDir.domain    = @(coor) isInUp(coor);
             sDir.direction = [2];
             sDir.value     = uVal*sin(angle);
             Dir3 = DirichletCondition(obj.mesh,sDir);             

             s.mesh = obj.mesh;
             s.dirichletFun = [Dir1 Dir2 Dir3];
             s.pointloadFun = [];
             s.periodicFun = [];
             obj.boundaryConditions = BoundaryConditions(s);
         end

        function createForceBendingConditions(obj,fVal)
            isInLeft = @(coord) (abs(coord(:,1) - min(coord(:,1)))< 1e-12);
            sDir.domain    = @(coor) isInLeft(coor);
            sDir.direction = [1,2];
            sDir.value     = 0;
            Dir1 = DirichletCondition(obj.mesh,sDir);

            isInRight = @(coord) (abs(coord(:,1) - max(coord(:,1)))< 1e-12);
            sNeum.domain    = @(coor) isInRight(coor);
            sNeum.direction = [2];
            sNeum.value     = fVal;
            Neum1 = DirichletCondition(obj.mesh,sNeum);

            s.mesh = obj.mesh;
            s.dirichletFun = [Dir1];
            s.pointloadFun = [Neum1];
            s.periodicFun = [];
            obj.boundaryConditions = BoundaryConditions(s);
        end

        function createDisplacementBendingConditions(obj,uVal)
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

         % function createLshapeDisplacementConditions(obj,uVal)
         %     isInDown = @(coor) (abs(coor(:,2) - min(coor(:,2)))  < 1e-12);
         %     sDir.domain    = @(coor) isInDown(coor);
         %     sDir.direction = [1,2];
         %     sDir.value     = 0;
         %     Dir1 = DirichletCondition(obj.mesh,sDir);
         % 
         %     isInTip = @(coor) (abs(coor(:,2)-(max(coor(:,2))+min(coor(:,2)))/2) < 1e-12) & (abs(coor(:,1)-max(coor(:,1))) < 30);
         %     sDir.domain    = @(coor) isInTip(coor);
         %     sDir.direction = [1];
         %     sDir.value     = 0;
         %     Dir2 = DirichletCondition(obj.mesh,sDir);
         % 
         %     sDir.domain    = @(coor) isInTip(coor);
         %     sDir.direction = [2];
         %     sDir.value     = uVal;
         %     Dir3 = DirichletCondition(obj.mesh,sDir);
         % 
         %     % Merge
         %     s.mesh = obj.mesh;
         %     s.dirichletFun = [Dir1 Dir2 Dir3];
         %     s.pointloadFun = [];
         %     s.periodicFun = [];
         %     obj.boundaryConditions = BoundaryConditions(s);
         % end
         % 
         % function createFiberMatrixDisplacementConditions(obj,uVal)
         %     % nodes = 1:obj.mesh.nnodes;
         %     % ndim = 2;
         %     % 
         %     % % Enforce fixed Dirichlet conditions to the down nodes
         %     % downSide = min(obj.mesh.coord(:,2));
         %     % isInDown = abs(obj.mesh.coord(:,2)-downSide) < 1e-12;
         %     % dirichletDown = zeros(ndim*length(nodes(isInDown)),3);
         %     % for i=1:ndim
         %     %     dirichletDown(i:2:end,1) = nodes(isInDown);
         %     %     dirichletDown(i:2:end,2) = i;
         %     % end
         %     % 
         %     % % Enforce fixed Dirichlet conditions to the fiber nodes
         %     % center = [(min(obj.mesh.coord(:,1))+max(obj.mesh.coord(:,1)))/2;
         %     %           (min(obj.mesh.coord(:,2))+max(obj.mesh.coord(:,2)))/2];
         %     % radius = 0.2;
         %     % isInCircle = ((obj.mesh.coord(:,1)-center(1)).^2 + (obj.mesh.coord(:,2)-center(2)).^2) ...
         %     %              < (radius^2 + 1e-5);
         %     % dirichletCircle = zeros(ndim*length(nodes(isInCircle)),3);
         %     % for i=1:ndim
         %     %     dirichletCircle(i:2:end,1) = nodes(isInCircle);
         %     %     dirichletCircle(i:2:end,2) = i;
         %     % end
         %     % 
         %     % % Enforce roller Dirichlet conditions to the top nodes
         %     % upSide  = max(obj.mesh.coord(:,2));
         %     % isInUp = abs(obj.mesh.coord(:,2)-upSide)< 1e-12;
         %     % dirichletUp   = zeros(ndim*length(nodes(isInUp)),3);
         %     % for i=1:ndim
         %     %     dirichletUp(i:2:end,1) = nodes(isInUp);
         %     %     dirichletUp(i:2:end,2) = i;
         %     % end
         %     % 
         %     % % Enforce displacement at the top
         %     % dirichletUp(2:2:end,3) = uVal;
         %     % 
         %     % % Merge
         %     % bc.dirichlet = [dirichletDown; dirichletUp; dirichletCircle];
         %     % bc.pointload = [];
         %     % obj.boundaryConditions = bc;
         % end

    end

end