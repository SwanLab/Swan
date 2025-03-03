classdef Perimeter3D < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        fileName
        topOptSet
        topOptProblem
    end
    
    methods (Access = public)
        
        function obj = Perimeter3D()
            obj.init();              
            obj.createTopOptSettings();            
            obj.solveProblem();
        end        
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.fileName = 'CantileverTetraPerimeterTotal';
           % obj.fileName = 'CantileverTetra';
        end
                
        function createTopOptSettings(obj)
          t = SettingsTopOptProblem(obj.fileName);            

          ls = obj.createLevelSet();         
          t.designVarSettings.initialCase = 'given';
          t.designVarSettings.creatorSettings.value = ls;               
          vF = t.incrementalSchemeSettings.targetParamsSettings.VfracFinal;            
          t.incrementalSchemeSettings.targetParamsSettings.VfracInitial = vF;     
          
          
          
          obj.topOptSet = t;              
        end
        
        function solveProblem(obj)
            obj.topOptProblem = TopOpt_Problem(obj.topOptSet);
            obj.topOptProblem.computeVariables();
            obj.topOptProblem.postProcess();
        end        
        
        function ls = createLevelSet(obj)
            ls = obj.loadLevelSet();
            %ls = obj.readLevelSet();
        end
        
        function ls = loadLevelSet(obj)
            %a = load('Output/CantileverTetra/DesignVariable320');
            a = load('/home/alex/Desktop/CantileverTetra/DesignVariable320');
            %a = load('Topology Optimization/Applications/PerimeterExperiments/LevelSetPerimeter3D');
            ls = a.x;
            m  = a.mesh;
            obj.plotByUnfittedMesh(m,ls);
         %   obj.plotSurfaceByInterpolation(m,ls);
        end
        
        function ls = readLevelSet(obj)
            iter = 500;
            fCase = 'CantileverTetraPerimeter';
            folder = '/home/alex/Desktop/CantileverTetraPerimeter/';            
            s.fileName = [fCase,num2str(iter)];
            s.folderPath = fullfile(folder);            
            wM = WrapperMshResFiles(s);
            wM.compute();
            ls = wM.dataRes.DesignVar1;
            obj.plotByUnfittedMesh(wM.mesh,ls);
            obj.plotSurfaceByInterpolation(wM.mesh,ls);
        end
        
        function plotSurfaceByInterpolation(obj,mesh,ls)
            x = mesh.coord(:,1);
            y = mesh.coord(:,2);
            z = mesh.coord(:,3);
            F = scatteredInterpolant(x,y,z,ls);
            xRange = linspace(min(x),max(x),75);
            yRange = linspace(min(y),max(y),75);
            zRange = linspace(min(z),max(z),75);
            [x,y,z] = meshgrid(xRange,yRange,zRange);
            ls = F(x,y,z);
            figure
            p = patch(isosurface(x,y,z,ls,0));
            isonormals(x,y,z,ls,p);
            set(p,'FaceColor','red','EdgeColor','none','FaceAlpha',1);
            daspect([1 1 1])
            view(3)
            axis tight
            grid on
            camlight; lighting phong
        end
        
        function plotByUnfittedMesh(obj,mesh,ls)
            s.backgroundMesh  = mesh;
            sB.backgroundMesh = mesh;
            sB.dimension = 1:3;
            sB.type = 'FromReactangularBox';
            bM = BoundaryMeshCreator.create(sB);
            s.boundaryMesh = bM.create();
            u = UnfittedMesh(s);
            u.compute(ls)
            figure
            u.plot();
        end
        
    end
    
end