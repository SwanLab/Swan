classdef ShFunc_NonSelfAdjoint_Compliance < Shape_Functional
    properties
        h_C_0; %compliance incial
        physicalProblem
        forces_adjoint
        adjointProblem
    end
    methods
        function obj = ShFunc_NonSelfAdjoint_Compliance(settings)
            obj@Shape_Functional(settings);
            switch settings.ptype
                case 'MACRO'
                    obj.physicalProblem = Physical_Problem(settings.filename);
                case 'MICRO'
                    obj.physicalProblem = Physical_Problem_Micro(settings.filename);
            end
            obj.physicalProblem.preProcess;
            
            obj.forces_adjoint = Preprocess.getBC_adjoint(settings.filename);
            obj.adjointProblem = Physical_Problem(settings.filename);
            
            [neumann_adj_dof,nuemann_adj_values] = obj.adjointProblem.dof.get_dof_conditions(obj.forces_adjoint,obj.adjointProblem.dim.nunkn);
            obj.adjointProblem.dof.neumann = neumann_adj_dof;
            obj.adjointProblem.dof.neumann_values = -nuemann_adj_values;
            
            obj.adjointProblem.preProcess;
        end
        function computef(obj,x,interpolation)
            rho = obj.filter.getP0fromP1(x);
            matProps = interpolation.computeMatProp(rho);
            
            %compute compliance
            obj.physicalProblem.setMatProps(matProps);
            obj.physicalProblem.computeVariables;
            
            strain = obj.physicalProblem.variables.strain;
            
            obj.adjointProblem.setMatProps(matProps);
            obj.adjointProblem.computeVariables;
            
            strain_adjoint = obj.adjointProblem.variables.strain;
            
            %compliance = d_u'*(-obj.adjointProblem.RHS);
            compliance = obj.physicalProblem.variables.d_u'*(-obj.adjointProblem.variables.fext);
            
            %compute gradient
            gradient_compliance = zeros(obj.physicalProblem.mesh.nelem,obj.physicalProblem.geometry.ngaus);
            for igaus = 1:obj.physicalProblem.geometry.ngaus
                for istre = 1:obj.physicalProblem.dim.nstre
                    for jstre = 1:obj.physicalProblem.dim.nstre
                        gradient_compliance(:,igaus) = gradient_compliance(:,igaus)+(squeeze(strain(igaus,istre,:)).*squeeze(matProps.dC(istre,jstre,:)).*squeeze(strain_adjoint(igaus,jstre,:)));
                    end
                end
            end
            
            %% !! NOTE: INVERSE ORDER THAN NORMAL COMPLIANCE, MISTAKE? !!
            
            if isempty(obj.h_C_0)
                obj.h_C_0 = compliance;
            else
                compliance = compliance/abs(obj.h_C_0);
                gradient_compliance = gradient_compliance/abs(obj.h_C_0);
            end
            
            gradient_compliance = obj.filter.getP1fromP0(gradient_compliance);
            gradient_compliance = obj.filter.Msmooth*gradient_compliance;
            
            obj.value = compliance;
            obj.gradient = gradient_compliance;
        end
    end
end
