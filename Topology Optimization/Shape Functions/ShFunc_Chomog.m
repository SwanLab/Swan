classdef ShFunc_Chomog < Shape_Functional
   
    properties (Access = public)
        Chomog
        tstress
        tstrain
        Chomog_Derivatives
        physicalProblem
        interpolation
        rho
        matProps
        Ch_star
    end
    
    methods (Access = public)
        
        function obj=ShFunc_Chomog(cParams)
            cParams.filterParams.quadratureOrder = 'LINEAR';            
            obj.init(cParams);
            obj.physicalProblem = FEM.create(cParams.filename);
            obj.physicalProblem.preProcess;
            obj.interpolation = Material_Interpolation.create(cParams.materialInterpolationParams);
        end
        
        function f = getPhysicalProblems(obj)
            f{1} = obj.physicalProblem;
        end
    end
    
    methods (Access = protected)
        function compute_Chomog_Derivatives(obj,x)
            obj.rho=obj.filter.getP0fromP1(x);
            obj.matProps=obj.interpolation.computeMatProp(obj.rho);
            
            nstre = obj.physicalProblem.element.getNstre();
            ngaus = obj.elemGradientSize.ngaus;
            nelem = obj.elemGradientSize.nelem;
            
            obj.Chomog_Derivatives = zeros(nstre,nstre,ngaus,nelem);
            for istreChomog = 1:nstre
                for jstreChomog = 1:nstre
                    for igaus=1:ngaus
                        for istre=1:nstre
                            for jstre = 1:nstre
                                obj.Chomog_Derivatives(istreChomog,jstreChomog,igaus,:) = ...
                                    squeeze(obj.Chomog_Derivatives(istreChomog,jstreChomog,igaus,:)) + ...
                                    (squeeze(obj.tstrain(istreChomog,igaus,istre,:))...
                                    .*squeeze(obj.matProps.dC(istre,jstre,:))...
                                    .*squeeze(obj.tstrain(jstreChomog,igaus,jstre,:)));
                            end
                        end
                    end
                end
            end
            
        end
        
        function r = derivative_projection_Chomog(obj,inv_matCh,alpha,beta)
            nstre = obj.physicalProblem.element.getNstre();
            ngaus = size(obj.tstrain,2);
            nelem = obj.physicalProblem.element.nelem;
            
            weights = alpha*beta';
            weights_inv = inv_matCh*weights*inv_matCh;
            DtC1 = zeros(ngaus,nelem);
            DtC = zeros(ngaus,nelem);
            for igaus=1:ngaus
                for i=1:nstre
                    for j=1:nstre
                        DtC1(igaus,:) = squeeze(obj.Chomog_Derivatives(i,j,igaus,:));
                        DtC(igaus,:) = DtC(igaus,:)- weights_inv(i,j)*DtC1(igaus,:);
                    end
                end
            end
            r = DtC;
        end
        
        function computePhysicalData(obj,x)
            obj.rho=obj.filter.getP0fromP1(x);
            obj.matProps=obj.interpolation.computeMatProp(obj.rho);
            obj.physicalProblem.setMatProps(obj.matProps);
            obj.physicalProblem.computeChomog;
            obj.Chomog = obj.physicalProblem.variables.Chomog;
            obj.tstrain = obj.physicalProblem.variables.tstrain;
            obj.tstress = obj.physicalProblem.variables.tstress;
        end
        
        function compute_Ch_star(obj,TOL,C_Cstar_case)
            E_plus = TOL.E_plus;
            E_minus = TOL.E_minus;
            nu_plus = TOL.nu_plus;
            nu_minus = TOL.nu_minus;
            
            switch C_Cstar_case
                case 'negative_poisson'
                    kappa_f = @(E,nu) E/2*(1-nu);
                    mu_f = @(E,nu) E/2*(1-nu);
                    
                    k_plus = kappa_f(E_plus,nu_plus);
                    mu_plus = mu_f(E_plus,nu_plus);
                    
                    k_minus = kappa_f(E_minus,nu_minus);
                    mu_minus = mu_f(E_minus,nu_minus);
                    
                    
                    nu = @(k,mu) (k-mu)/(k+mu);
                    E = @(k,mu) (4*k*mu)/(k+mu);
                    C = @(E,nu) (E/(1-nu*nu)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2]);
                    
                    kappa_nu_min = k_minus;
                    mu_nu_min = mu_plus;
                    
                    nu_min = nu(kappa_nu_min,mu_nu_min);
                    E_nu_min = E(kappa_nu_min,mu_nu_min);
                    C_nu_min = C(E_nu_min,nu_min);
                    obj.Ch_star = C_nu_min;
                    
                case 'nu_0_6' %From Sigmund Thesis% rho = 0.38
                    C = @(E,nu) (E/(1-nu*nu)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2]);
                    nu = -0.6;
                    E = (1-nu*nu)*0.04;
                    obj.Ch_star = C(E,nu);
                    
                case 'Seba' % Es=0.08; nus=-0.25
                    obj.Ch_star = [0.0853    -0.0213       0;
                        -0.0213    0.0853       0;
                        0         0    0.0533];
                    
                    
                case 'nu_0_8' %From Sigmund Thesis% rho = 0.25
                    C = @(E,nu) (E/(1-nu*nu)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2]);
                    nu = -0.8;
                    E = (1-nu*nu)*0.02;
                    obj.Ch_star = C(E,nu);
                    
                case 'Vfrac07'
                    obj.Ch_star =[
                        0.4256    0.2837         0
                        0.2837    0.4256         0
                        0         0    0.1419];
                    
                case 'Vfrac06'
                    obj.Ch_star =[
                        0.2909    0.1940         0
                        0.1940    0.2909         0
                        0         0    0.0970];
                    
                case 'Vfrac05'
                    obj.Ch_star =[
                        0.1892    0.1261         0
                        0.1261    0.1892         0
                        0         0    0.0631];
                    
                case 'Vfrac04'
                    obj.Ch_star =[
                        0.1141    0.0761         0
                        0.0761    0.1141         0
                        0         0    0.0380];
                    
                case 'Vfrac03'
                    obj.Ch_star =[
                        0.0611    0.0407         0
                        0.0407    0.0611         0
                        0         0    0.0204];
                    
                case 'Composite'
                    obj.Ch_star =[1 0.15 0;
                        0.15 0.5 0;
                        0 0 0.2];
                case 'HoneyComb'
                    obj.Ch_star =0.094*[1 0.75 0
                        0.75 1 0
                        0 0 0.125];
                case 'InvertedHoneyComb'
                    obj.Ch_star =0.08*[1 -0.5 0
                        -0.5 1 0
                        0 0 0.06];
                case 'AcousticZeroShearA'
                    obj.Ch_star =[1 1 0;
                        1 1 0;
                        0 0 0];
                case 'NegativePoiss06'
                    obj.Ch_star =0.04*[1 -0.6 0;
                        -0.6 1 0;
                        0 0 0.8];
            end
        end
    end
    
    methods (Access = protected, Static)
        function r = projection_Chomog(inv_matCh,alpha,beta)
            weights = alpha*beta';
            r = sum(sum(weights.*inv_matCh));
        end
    end
end
