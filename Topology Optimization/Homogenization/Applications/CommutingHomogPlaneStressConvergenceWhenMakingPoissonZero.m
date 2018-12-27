classdef CommutingHomogPlaneStressConvergenceWhenMakingPoissonZero < handle
    
    properties
        theta
        direction
        nu
        nuS
        weakTensor
        stiffTensor
        vphTensor
        vhpTensor
        errors
    end
    
    methods (Access = public)
        
        function obj = CommutingHomogPlaneStressConvergenceWhenMakingPoissonZero
            obj.init()
            obj.computeErrors()
            obj.plotErrors()                        
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.theta = 0.3;
            angle = pi/6;
            dir = [cos(angle) sin(angle) 0];
            obj.direction = Vector3D;
            obj.direction.setValue(dir);
        end
        
        function computeErrors(obj)
            n = 20;            
            %nus = (1/3).^[1:n-1];
            %obj.nuS = [nus,0];          
            obj.nuS = linspace(0,1/3,n);
            for iNu = 1:n
                obj.nu = obj.nuS(iNu);
                obj.computeWeakStiffTensor()
                obj.computeVoigtPlaneStressHomog();
                obj.computeVoigtHomogPlaneStressTensor
                obj.errors(iNu) = obj.computeError();
            end            
        end
        
        function computeWeakStiffTensor(obj)
            E1  = 1;
            nu1 = obj.nu;
            E0  = 1e-3*E1;
            nu0 = nu1;
            obj.stiffTensor = IsotropicConstitutiveTensor(E1,nu1);
            obj.weakTensor  = IsotropicConstitutiveTensor(E0,nu0);
        end
        
              
        function computeVoigtPlaneStressHomog(obj)
            c0     = obj.weakTensor;
            c1     = obj.stiffTensor;
            dir{1} = obj.direction;
            m1     = 1;
            seqHomog = VoigtPlaneStressHomogHomogenizer(c0,c1,dir,m1,obj.theta);
            obj.vphTensor  = seqHomog.getPlaneStressHomogenizedTensor();
        end
        
        function computeVoigtHomogPlaneStressTensor(obj)
            c0     = obj.weakTensor;
            c1     = obj.stiffTensor;
            dir{1} = obj.direction;
            m1     = 1;
            seqHomog = VoigtHomogPlaneStressHomogenizer(c0,c1,dir,m1,obj.theta);  
            obj.vhpTensor  = seqHomog.getPlaneStressHomogenizedTensor();
        end 
        
        function error = computeError(obj)
            c1 = obj.vhpTensor.getValue();
            c2 = obj.vphTensor.getValue();
            error = norm(c2-c1)/norm(c1);
        end
        
        function plotErrors(obj)
            figureID = figure;
            h = plot(obj.nuS,obj.errors,'-+');
            ca = get(figureID,'CurrentAxes');
            xl = get(ca,'xlabel');
            yl = get(ca,'ylabel');
            set(gca,'fontsize',20)
            set(xl,'string','\nu','FontSize',40)
            set(yl,'string','Error','FontSize',40)
            set(h,'LineWidth',5);
            plot_x0=1920;
            plot_y0=500;
            plot_width=600;
            plot_height=600;  
            set(figureID,'units','points','position',[plot_x0,plot_y0,plot_width,plot_height])
            output_figure_format = '-dpng';%'-dpdf'; %'-depsc'
            path = 'Topology Optimization/Homogenization/Applications/CommutingConvergence';
            print(figureID,path,output_figure_format);            
        end
    end
    
end

