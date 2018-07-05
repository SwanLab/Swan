classdef Objective_Function_AugLag < Objective_Function
    properties
        value_initial
        lambda
        penalty
        %         inactive_constr
    end
    methods
        function obj=Objective_Function_AugLag(settings)
            %             obj.lambda=zeros(1,settings.nconstr);
            %             obj.penalty=ones(1,settings.nconstr);
            
            % !! PATCH !!
            load(fullfile(pwd,'Allaire_ShapeOpt','meshSize'));
            load(fullfile(pwd,'Allaire_ShapeOpt','constraints_weights'));
            
            % !! Allaire's uses constant lambdas !!
            obj.lambda=[0 0]*(dx*dy);
            %             obj.lambda=[lagV lagP]*(dx*dy);
            obj.penalty=zeros(1,settings.nconstr);
        end
        function computeFunction(obj,cost,constraint)
%             obj.inactive_constr=-obj.lambda'>constraint.value;
%             constraint.value(obj.inactive_constr)=obj.lambda(obj.inactive_constr);
%             obj.value=cost.value + obj.lambda*constraint.value + 0.5*obj.penalty*(constraint.value.*constraint.value);
%             
%             !! PATCH !! --> Lambda's only used to compute value, NOT gradient
            load(fullfile(pwd,'Allaire_ShapeOpt','constraints_weights'));
            obj.value = cost.value + [lagV lagP]*constraint.value; % !! Used ShFunc_Perimeter instead of perimeter's Allaire function
        end
        function computeGradient(obj,cost,constraint)
            %             constraint.gradient(:,obj.inactive_constr)=0;
            obj.gradient=constraint.gradient*obj.lambda' + constraint.gradient*(obj.penalty'.*constraint.value) + cost.gradient;
        end
    end
end
