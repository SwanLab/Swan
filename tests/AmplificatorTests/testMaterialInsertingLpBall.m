classdef testMaterialInsertingLpBall < testShowingError
    
    properties (Access = protected)
        tol = 5e-2;
        testName = 'testMaterialInsertingLpBall';
    end
    
    properties (Access = private)
        density
        CtensorVademecum
        CtensorSIMPALL
        fileName
        vadVariables
        interpolator
        designVar
        plotter
    end
    
    
    methods (Access = public)
        
        function obj = testMaterialInsertingLpBall()
            obj.init();
            obj.computeConstitutiveFromVademecum();
            obj.computeDensityFromVademecum();
            obj.computeConstitutiveFromDensity();
            obj.plotConstitutiveTensors();
        end
        
    end
    
    methods (Access = protected)
        
        function computeError(obj)
            obj.error = 0;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.fileName = 'VademecumSmoothCorner';
            obj.createDesignVariableFromRandMxMy();
            obj.loadVademecumVariables();
            obj.createInterpolator();
        end
        
        function createDesignVariableFromRandMxMy(obj)
            a = 0.01;
            b = 0.99;
            obj.designVar = (b-a).*rand(2000,1) + a;
        end
        
        function loadVademecumVariables(obj)
            matFile   = [obj.fileName,'.mat'];
            file2load = fullfile('Output',obj.fileName,matFile);
            v = load(file2load);
            obj.vadVariables = v.d;
        end
        
        function createInterpolator(obj)
            sM.x = obj.vadVariables.domVariables.mxV;
            sM.y = obj.vadVariables.domVariables.myV;
            sI.mesh = StructuredMesh(sM);
            obj.interpolator = Interpolator(sI);
        end
        
        function computeConstitutiveFromVademecum(obj)
            s.vadVariables = obj.vadVariables;
            s.interpolator = obj.interpolator;
            ct = ConstitutiveTensorFromVademecum(s);
            obj.CtensorVademecum = ct.computeCtensor(obj.designVar);
        end
        
        function computeDensityFromVademecum(obj)
            s.vadVariables = obj.vadVariables;
            s.interpolator = obj.interpolator;
            dt = DensityFromVademecum(s);
            obj.density = dt.computeDensity(obj.designVar);
        end
        
        function computeConstitutiveFromDensity(obj)
            matProp.rho_plus = 1;
            matProp.rho_minus = 0;
            matProp.E_plus = 1;
            matProp.E_minus = 1e-3;
            matProp.nu_plus = 1/3;
            matProp.nu_minus = 1/3;
            d.constitutiveProperties = matProp;
            d.interpolation = 'SIMPALL';
            d.dim = '2D';
            d.typeOfMaterial = 'ISOTROPIC';
            mI = Material_Interpolation.create(d);
            material = mI.computeMatProp(obj.density);
            me = Material_Elastic_ISO_2D(size(obj.density,1));
            me.setProps(material);
            obj.CtensorSIMPALL = me.C;
        end
        
        function plotConstitutiveTensors(obj)
            obj.plotter.path = '/home/alex/Dropbox/Amplificators/Images/FourthOrderAmplificator/';
            obj.plotter.indeces = [1 1; 1 2; 1 3; 2 2; 2 3; 3 3];
            for index = 1:length(obj.plotter.indeces)
                obj.plotter.fig = figure(index);
                obj.plotter.istre = obj.plotter.indeces(index,1);
                obj.plotter.jstre = obj.plotter.indeces(index,2);
                obj.plotVademecumTensor();
                obj.plotSimpAllTensor();
                obj.addLegend();
                obj.printPlot();
            end
        end
        
        function plotSimpAllTensor(obj)
            C = obj.CtensorSIMPALL;
            c(:,1) = C(obj.plotter.istre,obj.plotter.jstre,:);
            obj.plotter.h{1} = plot(obj.density,c,'+');
        end
        
        function plotVademecumTensor(obj)
            C = obj.CtensorVademecum;
            c(:,1) = C(obj.plotter.istre,obj.plotter.jstre,:);
            obj.plotter.h{2} = plot(obj.density,c,'+');
            hold on
        end
        
        function addLegend(obj)
            obj.plotter.cleg = ['C_{',num2str(obj.plotter.istre),num2str(obj.plotter.jstre),'}'];
            leg1 = [obj.plotter.cleg,' Vademecum'];
            leg2 = [obj.plotter.cleg,' SIMPALL'];
            legend({leg1,leg2})
        end
        
        function printPlot(obj)
            plotName = [obj.plotter.cleg,'Comparison'];
            p = plotPrinter(obj.plotter.fig,obj.plotter.h);
            p.print([obj.plotter.path,plotName]);            
        end        
        
    end    
    
end