classdef ConnectivityComputer < handle

    properties (Access = public)

    end

    properties (Access = private)
        mesh
        designVariable
        conductivityInterpolator 
        massInterpolator
        filter
        filterAdjoint
        monitoringEigenModes
    end

    properties (Access = private)
        
    end

    methods (Access = public)

        function obj = ConnectivityComputer()
%             lambda1 = []
%             h = 1/100;
% %             for radius = 0:h:0.6 %0.5
%             radius_all = 0.0:h:0.6
%             for radius = radius_all %0.5
%                 radius
            obj.init();
            obj.createMesh();
            obj.createFilter();
            obj.createDesignVariable();  
            obj.createConductivityInterpolator();
            obj.createMassInterpolator();
            lambda = obj.computeEigenValueFunctional();
%                 lambda1 = [lambda1, lambda]
%             end
%             figure
%             semilogy(radius_all,lambda1,'-')
%             xlabel('Radius')
%             ylabel('First Eigenvalue')
%             grid on
% %             save('PDEFP.mat','lambda1')
%             save('LUMP.mat','lambda1')
% %             save('PDE.mat','lambda1')
        end

    end

    methods (Access = private)

        function init(obj)
            close all
        end

        function createMesh(obj)
            x1 = linspace(0,1,100);
            x2 = linspace(0,1,100);
            [xv,yv] = meshgrid(x1,x2);
            [F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            m = Mesh.create(s);
            obj.mesh = m;
        end

        function createDesignVariable(obj,radius) %,radius
%             s.type        = 'CircleInclusion';
%             s.xCoorCenter       = 0.5;
%             s.yCoorCenter       = 0.5;
%             s.radius      = 0.1;%  0.1; %
% % % % % 
%             s.type        = 'RingWithHorizontalCrack'; %'RingSDF';
%             s.innerRadius = 0.1;
%             s.outerRadius = 2.0;
%             s.xCoorCenter = 0.5;
%             s.yCoorCenter = 0.5;
%             s.crackWidth = 0.02;

            s.type = 'FourPerpendicularBars';
            s.leftBar_xMax = 0.35;   % right edge of left bar
            s.barWidth = 0.1;

            s.rightBar_xMin = 1 - s.leftBar_xMax;  % left edge of right bar
            s.bottomBar_yMax = s.leftBar_xMax ; % top edge of bottom bar
            s.topBar_yMin = s.rightBar_xMin;    % bottom edge of top bar

%             s.type = 'FourPerpendicularBarsWithCrack';
%             s.leftBar_xMax = 0.35;   % right edge of left bar
%             s.barWidth = 0.05;
% 
%             s.rightBar_xMin = 1 - s.leftBar_xMax;  % left edge of right bar
%             s.bottomBar_yMax = s.leftBar_xMax ; % top edge of bottom bar
%             s.topBar_yMin = s.rightBar_xMin;    % bottom edge of top bar
%             s.hCrack = 0.02;

                % % % % 
%             s.type        = 'RingSDF';
%             s.innerRadius = 0.1;
%             s.outerRadius = 0.15;
%             s.xCoorCenter = 0.5;
%             s.yCoorCenter = 0.5;

%             s.type        = 'ThreeRectanglesInclusion';
%             s.xSide1       = 0.3;
%             s.ySide1       = 0.3;
%             s.xCoorCenter1 = 0.35;
%             s.yCoorCenter1 = 0.5;
%             s.xSide2       = 0.1;
%             s.ySide2       = 0.1;
%             s.xCoorCenter2 = 0.6;
%             s.yCoorCenter2 = 0.5;
%             s.xSide3       = 3.0;
%             s.ySide3       = 0.1;
%             s.xCoorCenter3 = 0.5;
%             s.yCoorCenter3 = 0.2; 

            g              = GeometricalFunction(s);
            lsFun          = g.computeLevelSetFunction(obj.mesh);
            s.fun  = lsFun;
            s.mesh = obj.mesh;
            s.type = 'LevelSet';
            s.plotting = true;
            ls     = DesignVariable.create(s);
            obj.designVariable = ls;
        end

        function createFilter(obj)
%             s.filterType = 'PDE';
%             s.mesh  = obj.mesh;
%             s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
%             f = Filter.create(s);
%             f.updateEpsilon(1*obj.mesh.computeMinCellSize())
%             obj.filter = f;

            s.filterType = 'FilterAndProject';
%             s.filterType = 'CloseOperator'; %'FilterAndProject';
            s.mesh       = obj.mesh;
            s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
            s.filterStep = 'PDE';
            s.beta       = 100.0; % 1.0;
            s.eta        = 0.5;
            obj.filter = Filter.create(s);
            obj.filter.updateEpsilon(1*obj.mesh.computeMinCellSize())
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

        function [lambda] = computeEigenValueFunctional(obj)
            s.mesh = obj.mesh;
            s.designVariable = obj.designVariable;
            s.filter = obj.filter;
            s.boundaryConditions = obj.createEigenvalueBoundaryConditions();
            s.conductivityInterpolator = obj.conductivityInterpolator; 
            s.massInterpolator         = obj.massInterpolator; 
            s.dim = '2D';
            s.eigenModes = StiffnessEigenModesComputer(s);
            s.isCompl  = true;
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