classdef BcContinuumDamage < handle
    
    properties (Access = public)
        type
        bcValueSet
        valueSetLenght
        mesh
    end
    
    methods (Access = public)
        function obj = BcContinuumDamage (cParams)
            obj.init(cParams);
        end

        function bc = nextStep(obj,i)
            s.bcVal = obj.bcValueSet(i);
            bc = obj.bcSetType (s);
        end
    end

    methods (Access =  private)
        function init(obj,cParams) 
            obj.type = cParams.bcType;
            obj.bcValueSet = cParams.bcValueSet;
            obj.mesh = cParams.mesh;
            obj.valueSetLenght = size(obj.bcValueSet,2);
        end

        function  bc = bcSetType (obj, s)
            switch obj.type
                case 'displacementTraction'
                    isDown = @(coord) (abs(coord(:,2) - min(coord(:,2)))< 1e-12);
                    sDir.domain    = @(coor) isDown(coor);
                    sDir.direction = [1,2];
                    sDir.value     = 0;
                    Dir1 = DirichletCondition(obj.mesh,sDir);

                    isUp = @(coord) (abs(coord(:,2) - max(coord(:,2)))< 1e-12);
                    sDir.domain    = @(coor) isUp(coor);
                    sDir.direction = [2];
                    sDir.value     = s.bcVal;
                    Dir2 = DirichletCondition(obj.mesh,sDir);

                    s.mesh = obj.mesh;
                    s.dirichletFun = [Dir1 Dir2];
                    s.pointloadFun = [];
                    s.periodicFun = [];
                    bc = BoundaryConditions(s);

                case 'forceTraction'
                    isDown = @(coord) (abs(coord(:,2) - min(coord(:,2)))< 1e-12);
                    sDir.domain    = @(coor) isDown(coor);
                    sDir.direction = [1,2];
                    sDir.value     = 0;
                    Dir1 = DirichletCondition(obj.mesh,sDir);

                    isUp = @(coord) (abs(coord(:,2) - max(coord(:,2)))< 1e-12);
                    sNeum.domain    = @(coor) isUp(coor);
                    sNeum.direction = [2];
                    sNeum.value     = s.bcVal;
                    Neum1 = PointLoad(obj.mesh,sNeum);

                    s.mesh = obj.mesh;
                    s.dirichletFun = [Dir1];
                    s.pointloadFun = [Neum1];
                    s.periodicFun = [];
                    bc = BoundaryConditions(s);

                case 'displacementBending'
                    isLeft = @(coord) (abs(coord(:,1) - min(coord(:,1)))< 1e-12);
                    sDir.domain    = @(coor) isLeft(coor);
                    sDir.direction = [1,2];
                    sDir.value     = 0;
                    Dir1 = DirichletCondition(obj.mesh,sDir);

                    isRight = @(coord) (abs(coord(:,1) - max(coord(:,1)))< 1e-12);
                    sDir.domain    = @(coor) isRight(coor);
                    sDir.direction = [2];
                    sDir.value     = s.bcVal;
                    Dir2 = DirichletCondition(obj.mesh,sDir);

                    s.mesh = obj.mesh;
                    s.dirichletFun = [Dir1 Dir2];
                    s.pointloadFun = [];
                    s.periodicFun = [];
                    bc = BoundaryConditions(s);

                case 'forceBending'
                    isLeft = @(coord) (abs(coord(:,1) - min(coord(:,1)))< 1e-12);
                    sDir.domain    = @(coor) isLeft(coor);
                    sDir.direction = [1,2];
                    sDir.value     = 0;
                    Dir1 = DirichletCondition(obj.mesh,sDir);

                    isRight = @(coord) (abs(coord(:,1) - max(coord(:,1)))< 1e-12);
                    sNeum.domain    = @(coor) isRight(coor);
                    sNeum.direction = [2];
                    sNeum.value     = s.bcVal;
                    Neum1 = PointLoad(obj.mesh,sNeum);

                    s.mesh = obj.mesh;
                    s.dirichletFun = [Dir1];
                    s.pointloadFun = [Neum1];
                    s.periodicFun = [];
                    bc = BoundaryConditions(s);  
            end
        end
    end
end