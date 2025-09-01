classdef HyperelasticityBoundaryCreator < handle

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

        function obj = HyperelasticityBoundaryCreator(mesh,cParams)
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
            obj.type = cParams.type;
            obj.bcValues = cParams.bcValues;
            obj.step = 0;
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
            Neum1 = PointLoad(obj.mesh,sNeum);
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
            sDir.domain    = @(coor) isUp(coor);
            sDir.direction = [1];
            sDir.value     = 0;
            Dir2 = DirichletCondition(obj.mesh,sDir);

            sDir.domain    = @(coor) isUp(coor);
            sDir.direction = [2];
            sDir.value     = uVal;
            Dir3 = DirichletCondition(obj.mesh,sDir);

            s.mesh = obj.mesh;
            s.dirichletFun = [Dir1 Dir2 Dir3];
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

    end

end