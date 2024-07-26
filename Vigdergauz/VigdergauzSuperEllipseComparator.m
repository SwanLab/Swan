classdef VigdergauzSuperEllipseComparator < handle
    
    properties (Access = public)
        rhoDifferenceNorm
        frame
        mx
        my      
    end
    
    properties (Access = private)
        vigdergauzInputFile
        superEllipseInputFile
        settings
        mesh
        filter
        rhoDifference
        superellipseRho
        vigdergauzRho  
        quadrature
        dV
        qExponent
        vigdergauzLS
        superellipseLS
        superEllipseRatio
        volumeMicro         
    end
    
    methods (Access = public)
        
        function obj = VigdergauzSuperEllipseComparator()
            obj.init();
            obj.createFilter();
            obj.createQuadrature();
            obj.computeDV();   
        end
        
        function computeComparison(obj,txi,rho,q)
            obj.computeParameters(txi,rho,q)
            obj.createSuperEllipseLevelSet();
            obj.createVigdergauzLevelSet();
            obj.computeSuperEllipseDensity();
            obj.computeVigdergauzDensity();
            obj.computeDifference();
            obj.computeDifferenceNorm();            
        end
        
        function plotDifference(obj)
            obj.plotDensity(abs(obj.rhoDifference));
        end
                
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.vigdergauzInputFile = 'VigergauzLevelSetInput';
            obj.superEllipseInputFile = 'SuperEllipseLevelSetInput';
        end
        
        function createVigdergauzLevelSet(obj)
            obj.createSettings(obj.vigdergauzInputFile);
            obj.createMesh();
            obj.createVigdergauzSettings();
            obj.vigdergauzLS = obj.computeLevelSet();
        end
        
        function createSuperEllipseLevelSet(obj)
            obj.createSettings(obj.superEllipseInputFile);
            obj.createMesh();
            obj.createSuperEllipseSettings();
            obj.superellipseLS = obj.computeLevelSet();
        end
        
        function createSettings(obj,file)
            setting = Settings(file);
            translator = SettingsTranslator();
            translator.translate(setting);
            fileName = translator.fileName;
            obj.settings = SettingsTopOptProblem(fileName);
        end
        
        function vSet = createVigdergauzSettings(obj)
            vSet = obj.settings;
            vSet.designVarSettings.mesh = obj.mesh;
            lsS = vSet.designVarSettings.levelSetCreatorSettings;
            lsS.vigdergauzDataBase.superEllipseRatio = obj.superEllipseRatio;
            lsS.vigdergauzDataBase.volumeMicro       = obj.volumeMicro;
            vSet.designVarSettings.levelSetCreatorSettings = lsS;
            obj.settings = vSet;
        end
        
        function computeParameters(obj,txi,rho,q)
            obj.superEllipseRatio = tan(txi);
            obj.volumeMicro = rho;  
            se = SuperEllipseParamsRelator();
            obj.mx = se.mx(txi,rho,q);
            obj.my = se.my(txi,rho,q);
            obj.qExponent  = q;                        
        end
        
        function createSuperEllipseSettings(obj)
            set = obj.settings;
            set.designVarSettings.mesh = obj.mesh;
            lsS = set.designVarSettings.levelSetCreatorSettings;
            lsS.widthH  = obj.mx;
            lsS.widthV = obj.my;
            lsS.pnorm  = obj.qExponent;
            set.designVarSettings.levelSetCreatorSettings = lsS;
            obj.settings = set;
        end
        
        function createMesh(obj)
            s.coord  = obj.settings.designVarSettings.femData.coord;
            s.connec = obj.settings.designVarSettings.femData.connec;
            obj.mesh = Mesh_Total(s);
        end
        
        function ls = computeLevelSet(obj)
            s = obj.settings.designVarSettings;
            s.scalarProductSettings.epsilon = 1e-3;
            designV = DesignVariable.create(s);
            ls = designV.value;
        end
        
        function createFilter(obj)
            obj.createSettings(obj.superEllipseInputFile);
            obj.createMesh();            
            designV.mesh = obj.mesh;
            designV.type = 'LevelSet';
            filterParams = obj.settings.costSettings.shapeFuncSettings{1}.filterParams;
            filterParams.designVar = designV;
            filterParams.quadratureOrder = 'LINEAR';
            obj.filter = Filter.create(filterParams);
%             obj.filter.preProcess();
        end
        
        function createQuadrature(obj)
            q = Quadrature.set('TRIANGLE');
            q.computeQuadrature('LINEAR'); 
            obj.quadrature = q;
        end        
        
        function computeDV(obj)            
            q = obj.quadrature;
            obj.dV = obj.mesh.computeDvolume(q);
        end
        
        function computeSuperEllipseDensity(obj)
            rho = obj.filter.getP0fromP1(obj.superellipseLS);
            obj.superellipseRho = rho;
        end
        
        function computeVigdergauzDensity(obj)
            rho = obj.filter.getP0fromP1(obj.vigdergauzLS);
            obj.vigdergauzRho = rho;
        end
        
        function computeDifference(obj)
            dif = obj.superellipseRho - obj.vigdergauzRho;
            obj.rhoDifference = dif;    
        end
        
        function normF = computeNorm(obj,f)
            dv(:,1) = obj.dV;
            q = obj.quadrature;
            int = 0;            
            for igaus = 1:q.ngaus
               int = int + sum(f.*f.*dv);
            end    
            normF = int;
        end
        
        function computeDifferenceNorm(obj)
            rhoDiffNorm = obj.computeNorm(obj.rhoDifference);
            rhoSeNorm   = obj.computeNorm(obj.superellipseRho);           
            relNorm = rhoDiffNorm/rhoSeNorm; 
            obj.rhoDifferenceNorm = relNorm;
        end
        
        
        function plotDensity(obj,density)
            [x,y] = obj.computeXYgaussPoints();
            v = density;      
            xg = linspace(min(x),max(x),500);
            yg = linspace(min(y),max(y),500);
            [xq,yq] = meshgrid(xg,yg);
            vq = griddata(x,y,v,xq,yq);
            figureID = figure(1);
            h = surf(xq,yq,vq);            
            set(h,'LineStyle','none')
            view([0 0 1]);
            c = gray;
            c = flipud(c);
            colormap(c)
            vTitle = ['Volume = ',num2str(obj.rhoDifferenceNorm*100),'%'];
            qTitle = ['q = ',num2str(obj.qExponent)];
            fName = [vTitle, ' {and} ',qTitle];
            title(fName);
            drawnow
            c = contourPrinter(figureID);
            c.print([fName,'.png'])
            obj.frame = getframe(gcf);
        end
        
        function [x,y] = computeXYgaussPoints(obj)
            connec = obj.mesh.connec;
            coord  = obj.mesh.coord;
            nnode = size(connec,2);
            nelem = size(connec,1);
            q = obj.quadrature;
            x = zeros(nelem,q.ngaus);
            y = zeros(nelem,q.ngaus);            
            xpg = q.posgp;            
            inter = Interpolation.create(obj.mesh.type,'LINEAR');
            for igaus = 1:q.ngaus
                for i = 1:nnode
                    shape = inter.computeShapeFunctions(xpg(:,igaus));
                    xp = coord(connec(:,i),1);
                    yp = coord(connec(:,i),2);                    
                    x(:,igaus) = x(:,igaus) + shape(1,igaus)*xp;
                    y(:,igaus) = y(:,igaus) + shape(2,igaus)*yp;
                end
            end                   
        end        
        
    end
    
end