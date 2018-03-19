classdef ShFunc_Chomog_EnforceCh_CCstar < ShFunc_Chomog_EnforceCh
    properties
        component
        epsilon
        initial_value = 1;
    end
    methods
        function obj=ShFunc_Chomog_EnforceCh_CCstar(settings,n)
            obj@ShFunc_Chomog_EnforceCh(settings);
            obj.Ch_star = obj.compute_Ch_star(settings.TOL);
            obj.selectiveC_Cstar=settings.selectiveC_Cstar;
            obj.component = n;
            obj.epsilon = settings.epsilon_isotropy_final;
        end
        function computef(obj,x)
            obj.physicalProblem.computeChomog;
            obj.setPhysicalData(obj.physicalProblem.variables);
            %Cost
            Ch_star_div = obj.Ch_star;
            Ch_star_div (abs(Ch_star_div) < 1e-5) = 1;
            C_C = (obj.Chomog - obj.Ch_star)./Ch_star_div;
            
            % C-C*
            sq2 = sqrt(2);
            weights = [1,1,1,sq2,sq2,sq2]';
            costval = weights.*[C_C(1,1); ...
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
            
            gradientval = zeros(obj.physicalProblem.element.nelem,neq);
            obj.compute_Chomog_Derivatives(x);
            for i = 1:neq
                C_C = selectiveC_Cstar(:,:,i)./selective_Ch_star_div(i);                
                DtC1 = zeros(obj.physicalProblem.element.quadrature.ngaus,obj.physicalProblem.element.nelem);
                DtC = zeros(obj.physicalProblem.element.quadrature.ngaus,obj.physicalProblem.element.nelem);
                for igaus=1:obj.physicalProblem.element.quadrature.ngaus
                    for a=1:obj.physicalProblem.element.nstre
                        for b=1:obj.physicalProblem.element.nstre
                            DtC1(igaus,:) = squeeze(obj.Chomog_Derivatives(a,b,igaus,:));
                            DtC(igaus,:) = DtC(igaus,:) + C_C(a,b)*DtC1(igaus,:);
                        end
                    end
                end
                gradientval(:,i) = weights(i)*DtC;
            end
            obj.value = (costval(obj.component))^2 - obj.epsilon;
            obj.gradient = 2*costval(obj.component)*gradientval(:,obj.component);
            obj.value = obj.value/obj.initial_value;
            obj.gradient = obj.gradient/obj.initial_value;
            obj.passFilter();
        end
    end
    
    methods(Access = public)
        function obj =  setEpsilon(obj,r)
            obj.epsilon = r;
        end
    end
end