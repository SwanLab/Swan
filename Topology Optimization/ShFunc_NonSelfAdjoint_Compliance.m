classdef ShFunc_NonSelfAdjoint_Compliance < ShFunc_Compliance
    properties
        forces_adjoint
        adjointProblem
    end
    methods
        function obj = ShFunc_NonSelfAdjoint_Compliance(settings)
            obj@ShFunc_Compliance(settings);
            
            obj.forces_adjoint = Preprocess.getBC_adjoint(settings.filename);
            obj.adjointProblem = Physical_Problem(settings.filename);
            
            [neumann_adj_dof,nuemann_adj_values] = obj.adjointProblem.dof.get_dof_conditions(obj.forces_adjoint,obj.adjointProblem.dof.nunkn);
            obj.adjointProblem.dof.neumann = neumann_adj_dof;
            obj.adjointProblem.dof.neumann_values = -nuemann_adj_values;
            
            obj.adjointProblem.preProcess;
        end
        function computef(obj,x)
            obj.rho = obj.filter.getP0fromP1(x);
            obj.matProps = obj.interpolation.computeMatProp(obj.rho);
            obj.adjointProblem.setMatProps(obj.matProps);
            obj.adjointProblem.computeVariables;
            
            computef_CORE(obj);
        end
        
        function compliance = computeCompliance(obj)
            compliance = obj.physicalProblem.variables.d_u'*(-obj.adjointProblem.variables.fext);
        end
        
        function gradient_compliance = updateGradient(obj,igaus,istre,jstre)
            gradient_compliance = (squeeze(obj.physicalProblem.variables.strain(igaus,istre,:)).*squeeze(obj.matProps.dC(istre,jstre,:)).*squeeze(obj.adjointProblem.variables.strain(igaus,jstre,:)));
        end
    end
end
