classdef ShFunc_Chomog_EnforceCh< ShFunc_Chomog
    properties (Access = protected)
        Ch_star
        %selectiveC_Cstar
    end
    methods
        function obj=ShFunc_Chomog_EnforceCh(settings)
            obj@ShFunc_Chomog(settings);
        end
    end
    methods (Access = protected)
        function obj = passFilter(obj)
            mass=obj.filter.Msmooth;
            g = obj.gradient;
            obj.gradient = zeros(size(mass,1),size(g,2));
            for t=1:size(obj.gradient,2)
                gradient=obj.filter.getP1fromP0(g(:,t));
                obj.gradient(:,t) = mass*gradient;
            end
            %             if isempty(obj.h_C_0)
            %                 obj.h_C_0 = obj.value;
            %             end
            %             obj.value = obj.value/abs(obj.h_C_0);
            %             gradient=gradient/abs(obj.h_C_0);
            %             obj.h_C_0 = costfunc;
        end
        function computeCCstar(obj,x)
            %Cost
            Ch_star_div = obj.Ch_star;
            Ch_star_div (abs(Ch_star_div) < 1e-5) = 1;
            C_C = (obj.Chomog - obj.Ch_star)./Ch_star_div;
            
            % C-C*
            sq2 = sqrt(2);
            weights = [1,1,1,sq2,sq2,sq2]';
            obj.value = weights.*[C_C(1,1); ...
                C_C(2,2); ...
                C_C(3,3); ...
                C_C(2,3); ...
                C_C(1,3); ...
                C_C(1,2)];
            
            %Gradient
            
            neq = 6;
            selectiveC_Cstar = zeros(3,3,neq);
            
            % Eqn 1
            selectiveC_Cstar(1,1,1) = 1;
            selective_Ch_star_div(1) = Ch_star_div(1,1);
            
            % Eqn 2
            selectiveC_Cstar(2,2,2) = 1;
            selective_Ch_star_div(2) = Ch_star_div(2,2);
            
            % Eqn 3
            selectiveC_Cstar(3,3,3) = 1;
            selective_Ch_star_div(3) = Ch_star_div(3,3);
            
            % Eqn 4
            selectiveC_Cstar(2,3,4) = 1;
            selective_Ch_star_div(4) = Ch_star_div(2,3);
            
            % Eqn 5
            selectiveC_Cstar(1,3,5) = 1;
            selective_Ch_star_div(5) = Ch_star_div(1,3);
            
            % Eqn 6
            selectiveC_Cstar(1,2,6) = 1;
            selective_Ch_star_div(6) = Ch_star_div(1,2);
            
            %Gradient
            obj.gradient = zeros(obj.physicalProblem.mesh.nelem,neq);
            obj.compute_Chomog_Derivatives(x);
            for i = 1:neq
                C_C = selectiveC_Cstar(:,:,i)./selective_Ch_star_div(i);
                DtC1 = zeros(obj.physicalProblem.geometry.ngaus,obj.physicalProblem.mesh.nelem);
                DtC = zeros(obj.physicalProblem.geometry.ngaus,obj.physicalProblem.mesh.nelem);
                for igaus=1:obj.physicalProblem.geometry.ngaus
                    for a=1:obj.physicalProblem.dim.nstre
                        for b=1:obj.physicalProblem.dim.nstre
                            DtC1(igaus,:) = squeeze(obj.Chomog_Derivatives(a,b,igaus,:));
                            DtC(igaus,:) = DtC(igaus,:) + C_C(a,b)*DtC1(igaus,:);
                        end
                    end
                end
                obj.gradient(:,i) = weights(i)*DtC;
            end
        end
        function compute_Ch_star(obj,TOL)
            E_plus = TOL.E_plus;
            E_minus = TOL.E_minus;
            nu_plus = TOL.nu_plus;
            nu_minus = TOL.nu_minus;
            
            C_Cstar_case = 'Vfrac03';
            
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
                    
                case 'Vfrac03'
                    obj.Ch_star = [0.0611    0.0407         0
                                0.0407    0.0611         0
                                    0         0    0.0204];
            end
        end
    end
end
