classdef Bc_ContinuumDamage < handle
 
    properties (Access = public)
        type
        bcValueSet
        ValueSetLenght
        index
        mesh

        bc
    end
    
    properties (Access = private)

    end
    
    methods (Access = public)
        function obj = Bc_ContinuumDamage (cParams)
            obj.init(cParams);
        end

        function nextStep (obj,i)
            obj.index = i;
            s.bcVal = obj.bcValueSet(obj.index);
            obj.bc = obj.bcSetType (s);
        end
    
    end

    methods (Access =  private)

        function init(obj,cParams) 
            obj.type = cParams.bcType;
            obj.bcValueSet = cParams.bcValueSet;
            obj.mesh = cParams.mesh;
            obj.index = 1;
            obj.ValueSetLenght = size(obj.bcValueSet,2);

            s.bcVal = obj.bcValueSet(obj.index);
            obj.bc = obj.bcSetType (s);

        end

        function  bc = bcSetType (obj, s)
            switch obj.type
                case 'displacementTraction'
                    isInDown = @(coord) (abs(coord(:,2) - min(coord(:,2)))< 1e-12);
                    sDir.domain    = @(coor) isInDown(coor);
                    sDir.direction = [1,2];
                    sDir.value     = 0;
                    Dir1 = DirichletCondition(obj.mesh,sDir);

                    isInUp = @(coord) (abs(coord(:,2) - max(coord(:,2)))< 1e-12);
                    sDir.domain    = @(coor) isInUp(coor);
                    sDir.direction = [2];
                    sDir.value     = s.bcVal;
                    Dir2 = DirichletCondition(obj.mesh,sDir);

                    s.mesh = obj.mesh;
                    s.dirichletFun = [Dir1 Dir2];
                    s.pointloadFun = [];
                    s.periodicFun = [];
                    bc = BoundaryConditions(s);

                case 'forceTraction'
                    isInDown = @(coord) (abs(coord(:,2) - min(coord(:,2)))< 1e-12);
                    sDir.domain    = @(coor) isInDown(coor);
                    sDir.direction = [1,2];
                    sDir.value     = 0;
                    Dir1 = DirichletCondition(obj.mesh,sDir);

                    isInUp = @(coord) (abs(coord(:,2) - max(coord(:,2)))< 1e-12);
                    sNeum.domain    = @(coor) isInUp(coor);
                    sNeum.direction = [2];
                    sNeum.value     = s.bcVal;
                    Neum1 = PointLoad(obj.mesh,sNeum);


                    s.mesh = obj.mesh;
                    s.dirichletFun = [Dir1];
                    s.pointloadFun = [Neum1];
                    s.periodicFun = [];
                    bc = BoundaryConditions(s);

                case 'displacementBending'
                    isInLeft = @(coord) (abs(coord(:,1) - min(coord(:,1)))< 1e-12);
                    sDir.domain    = @(coor) isInLeft(coor);
                    sDir.direction = [1,2];
                    sDir.value     = 0;
                    Dir1 = DirichletCondition(obj.mesh,sDir);

                    isInRight = @(coord) (abs(coord(:,1) - max(coord(:,1)))< 1e-12);
                    sDir.domain    = @(coor) isInRight(coor);
                    sDir.direction = [2];
                    sDir.value     = s.bcVal;
                    Dir2 = DirichletCondition(obj.mesh,sDir);

                    s.mesh = obj.mesh;
                    s.dirichletFun = [Dir1 Dir2];
                    s.pointloadFun = [];
                    s.periodicFun = [];
                    bc = BoundaryConditions(s);

                case 'forceBending'
                    isInLeft = @(coord) (abs(coord(:,1) - min(coord(:,1)))< 1e-12);
                    sDir.domain    = @(coor) isInLeft(coor);
                    sDir.direction = [1,2];
                    sDir.value     = 0;
                    Dir1 = DirichletCondition(obj.mesh,sDir);

                    isInRight = @(coord) (abs(coord(:,1) - max(coord(:,1)))< 1e-12);
                    sDir.domain    = @(coor) isInRight(coor);
                    sDir.direction = [2];
                    sDir.value     = s.bcVal;
                    Dir2 = PointLoad(obj.mesh,sDir);

                    s.mesh = obj.mesh;
                    s.dirichletFun = [Dir1 Dir2];
                    s.pointloadFun = [];
                    s.periodicFun = [];
                    bc = BoundaryConditions(s);  
            end
        end



    end
end