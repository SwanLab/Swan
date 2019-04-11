classdef testMaterialInsertingLpBall < testShowingError
    
    properties (Access = protected)
        tol = 5e-2;
        testName = 'testMaterialInsertingLpBall';
    end
    
    properties (Access = private)
        density
        CtensorVademecum
        CtensorSIMPALL
        CtensorDif
    end
    
    
    methods (Access = public)
        
        function obj = testMaterialInsertingLpBall()
            obj.computeConstitutiveFromVademecum();
            obj.computeConstitutiveFromDensity();
            obj.computeConstitutiveTensorDifference();
        end
        
    end
    
    methods (Access = protected)
        
        function computeError(obj)
            obj.error = norm(obj.CtensorDif(:));
        end
        
    end
    
    methods (Access = private)
        
        function computeConstitutiveFromVademecum(obj)
            a = 0.01;
            b = 0.99;
            designVar = (b-a).*rand(5,2) + a;

            
            s = SettingsConstitutiveTensorFromVademecum();
            s.fileName = 'VademecumSmoothCorner';

            ct = ConstitutiveTensorFromVademecum(s);                       
            obj.CtensorVademecum = ct.computeCtensor(designVar);
            
            dt = DensityFromVademecum(s);
            obj.density = dt.computeDensity(designVar);
            
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
        
        function computeConstitutiveTensorDifference(obj)
            Cv = obj.CtensorVademecum;
            Cs = obj.CtensorSIMPALL;
            indeces = [1 1; 1 2; 1 3; 2 2; 2 3; 3 3];
            path = '/home/alex/Dropbox/Amplificators/Images/FourthOrderAmplificator/';            
            for index = 1:length(indeces)
                    istre = indeces(index,1);
                    jstre = indeces(index,2);
                    cv(:,1) = Cv(istre,jstre,:);
                    cs(:,1) = Cs(istre,jstre,:);
                    fig = figure();
                    h{1} = plot(obj.density,cv,'+');
                    hold on
                    h{2} = plot(obj.density,cs,'+');                    
                    cleg = ['C_{',num2str(istre),num2str(jstre),'}'];
                    leg1 = [cleg,' Vademecum'];
                    leg2 = [cleg,' SIMPALL'];
                    legend({leg1,leg2})
                    plotName = [cleg,'Comparison'];
                    p = plotPrinter(fig,h);
                    p.print([path,plotName]);
            end
            obj.CtensorDif = abs(Cv - Cs)./max(1,Cs);
        end
        
    end
    
    
end