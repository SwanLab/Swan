classdef ShFunc_NonSelfAdjoint_Compliance < ShFunc_Compliance
    properties
        forces_adjoint
        adjointProb
    end
    methods
        function obj = ShFunc_NonSelfAdjoint_Compliance(settings)
            obj@ShFunc_Compliance(settings);
            
            obj.forces_adjoint = Preprocess.getBC_adjoint(settings.filename);
            obj.adjointProb = FEM.create(settings.filename);
            
            [neumann_adj_dof,neumann_adj_values] = obj.adjointProb.dof.get_dof_conditions(obj.forces_adjoint,obj.adjointProb.dof.nunkn);
            obj.adjointProb.dof.neumann = neumann_adj_dof;
            obj.adjointProb.dof.neumann_values = -neumann_adj_values;
            
            obj.adjointProb.preProcess;
        end
        function computef(obj,x)
            obj.rho = obj.filter.getP0fromP1(x);
            obj.matProps = obj.interpolation.computeMatProp(obj.rho);
            obj.adjointProb.setMatProps(obj.matProps);
            obj.adjointProb.computeVariables;
            
            computef_CORE(obj);
        end
        
        function compliance = computeCompliance(obj)
            compliance = obj.physProb.variables.d_u'*(-obj.adjointProb.variables.fext);
        end
        
        function gradient_compliance = updateGradient(obj,igaus,istre,jstre)
            gradient_compliance = (squeeze(obj.physProb.variables.strain(igaus,istre,:)).*squeeze(obj.matProps.dC(istre,jstre,:)).*squeeze(obj.adjointProb.variables.strain(igaus,jstre,:)));
        end
    end
end
