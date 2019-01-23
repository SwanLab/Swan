classdef ShFunc_Chomog_EnforceCh< ShFunc_Chomog
    properties (Access = protected)
    end
    
    methods
        function obj=ShFunc_Chomog_EnforceCh(settings)
            obj.init(settings);
        end
    end
    
    methods (Access = protected)
        function obj = passFilter(obj)
            mass=obj.Msmooth;
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
            nelem = obj.physicalProblem.element.nelem;
            ngaus = size(obj.tstrain,2);
            nstre = obj.physicalProblem.element.nstre;
            
            obj.gradient = zeros(nelem,neq);
            obj.compute_Chomog_Derivatives(x);
            for i = 1:neq
                C_C = selectiveC_Cstar(:,:,i)./selective_Ch_star_div(i);
                DtC1 = zeros(ngaus,nelem);
                DtC = zeros(ngaus,nelem);
                for igaus=1:ngaus
                    for a=1:nstre
                        for b=1:nstre
                            DtC1(igaus,:) = squeeze(obj.Chomog_Derivatives(a,b,igaus,:));
                            DtC(igaus,:) = DtC(igaus,:) + C_C(a,b)*DtC1(igaus,:);
                        end
                    end
                end
                obj.gradient(:,i) = weights(i)*DtC;
            end
        end
    end
end