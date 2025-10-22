classdef ConnectivityComputer3D < handle

    properties (Access = public)

    end

    properties (Access = private)
        mesh
        levelSet
        characteristicFunction
        designVariable
        materialInterpolation
        filter
        filterConnect
        filterAdjointConnect
        monitoringEigenModes
    end

    properties (Access = private)
        
    end

    methods (Access = public)

        function obj = ConnectivityComputer3D()
            obj.init();
            obj.createMesh();
            obj.createFilter();
            obj.createDesignVariable();  
            obj.createConductivityInterpolator();
            obj.createMassInterpolator();
            obj.computeEigenValueFunctional();
        end

    end

    methods (Access = private)

        function init(obj)
            close all
        end

        function createMesh(obj)
            obj.mesh = HexaMesh(1.0,1.0,1.0,40,40,40); %20,20,20);
        end

        function createDesignVariable(obj)
%             s.type        = 'RectangleInclus10^-3}, ion';
%             s.xSide       = 0.5;
%             s.ySide       = 0.5;
%             s.xCoorCenter = 0.5;
%             s.yCoorCenter = 0.5;

            s.type        = 'ThreePrisms';
            s.xSide1       = 0.2;
            s.ySide1       = 0.2;
            s.zSide1       = 0.2;
            s.xCoorCenter1 = 0.3;
            s.yCoorCenter1 = 0.7;
            s.zCoorCenter1 = 0.7;
            s.xSide2       = 2.0;
            s.ySide2       = 0.2;
            s.zSide2       = 0.2;
            s.xCoorCenter2 = 0.5;
            s.yCoorCenter2 = 0.2;
            s.zCoorCenter2 = 0.2;
            s.xSide3       = 0.3;
            s.ySide3       = 0.3;
            s.zSide3       = 0.3;
            s.xCoorCenter3 = 0.6;
            s.yCoorCenter3 = 0.4; 
            s.zCoorCenter3 = 0.6; 
            g             = GeometricalFunction(s);
            lsFun          = g.computeLevelSetFunction(obj.mesh);
            s.fun  = lsFun;
            s.mesh = obj.mesh;
            s.type = 'LevelSet';
            s.plotting = true;
            ls     = DesignVariable.create(s);
            obj.designVariable = ls;
        end

        function createFilter(obj)
            s.filterType = 'LUMP';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
            obj.filter = f;
        end        

        function createConductivityInterpolator(obj) 
            s.interpolation  = 'SimpAllThermal';
            s.f0   = 1e-3;                                             
            s.f1   = 1;  
            s.dim  = '2D';
            a = MaterialInterpolator.create(s);
            obj.conductivityInterpolator = a;            
        end 

        function createMassInterpolator(obj)
            s.interpolation  = 'SIMPThermal';                              
            s.f0   = 1e-3;
            s.f1   = 1;
            s.pExp = 1;
            a = MaterialInterpolator.create(s);
            obj.massInterpolator = a;            
        end      

        function computeEigenValueFunctional(obj)
            s.mesh = obj.mesh;
            s.designVariable = obj.designVariable;
            s.filter = obj.filter;
            s.boundaryConditions = obj.createEigenvalueBoundaryConditions();
            s.eigenModes = StiffnessEigenModesComputer(s);
            mE = MinimumEigenValueFunctional(s);
            [lambda, dlambda] = mE.computeFunctionAndGradient(obj.designVariable);  
        end

        function  bc = createEigenvalueBoundaryConditions(obj)
            xMin    = min(obj.mesh.coord(:,1));
            yMin    = min(obj.mesh.coord(:,2));
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            isLeft  = @(coor) abs(coor(:,1))==xMin;
            isRight = @(coor) abs(coor(:,1))==xMax;
            isFront = @(coor) abs(coor(:,2))==yMin;
            isBack = @(coor) abs(coor(:,2))== yMax;
            isDir   = @(coor) isLeft(coor) | isRight(coor) | isFront(coor) | isBack(coor);  
            sDir{1}.domain    = @(coor) isDir(coor);
            sDir{1}.direction = 1;
            sDir{1}.value     = 0;
            sDir{1}.ndim = 1;
            
            dirichletFun = [];
            for i = 1:numel(sDir)
                dir = DirichletCondition(obj.mesh, sDir{i});
                dirichletFun = [dirichletFun, dir];
            end
            s.dirichletFun = dirichletFun;
            s.pointloadFun = [];

            s.periodicFun  = [];
            s.mesh         = obj.mesh;
            bc = BoundaryConditions(s);  
        end

    end

end