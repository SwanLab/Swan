classdef ShFunc_StressNorm < ShFunWithElasticPdes
    
    properties (Access = private)
        adjointProb
    end
    
    methods (Access = public)
        
        function obj = ShFunc_StressNorm(cParams)
            cParams.filterParams.quadratureOrder = 'LINEAR';            
            obj.init(cParams);     
            obj.createEquilibriumProblem(cParams.filename);
           % obj.createAdjointProblem(cParams.filename)            
        end
        
        function f = getPhysicalProblems(obj)
            f{1} = obj.physicalProblem;
            f{2} = obj.adjointProb;
        end
        
    end
    
    methods (Access = protected)

        function computeFunctionValue(obj)
            phy = obj.physicalProblem;
            u = phy.variables.d_u;
            f = phy.variables.fext;
            c2 = -f'*u;   
            
            
            dvolum  = phy.geometry.dvolu;
            strain = phy.variables.strain;
            C = phy.element.material.C;
            ngaus = phy.element.quadrature.ngaus;
            nstre = size(strain,2);
            c = 0;
            for igaus = 1:ngaus
                strainG = squeeze(strain(igaus,:,:));
                dV = dvolum(:,igaus);
                for istre = 1:nstre
                    sI(:,1) = strainG(istre,:);
                    for jstre = 1:nstre
                        sJ(:,1) = strainG(jstre,:);
                        Cij = squeeze(C(istre,jstre,:));
                        csum = sI.*Cij.*sJ.*dV;
                        c = c + sum(csum);
                    end
                end
            end            
            
            obj.value = c
            
            
            
            % Sum with stress + Amplificators
            %
            
        end
        
        function updateHomogenizedMaterialProperties(obj)
            obj.filterDesignVariable();
            obj.homogenizedVariablesComputer.computeCtensor(obj.regDesignVariable);
            obj.homogenizedVariablesComputer.computePtensor(obj.regDesignVariable);
        end
                
        
        function g = updateGradient(obj,igaus,istre,jstre)
            eu   = obj.physicalProblem.variables.strain;
            ev   = obj.adjointProb.variables.strain;
            eu_i = squeeze(eu(igaus,istre,:));
            ev_j = squeeze(ev(igaus,jstre,:)); 
            g = zeros(length(eu_i),1,obj.nVariables);            
            for ivar = 1:obj.nVariables
                dCij = squeeze(obj.homogenizedVariablesComputer.dC(istre,jstre,ivar,:));            
                g(:,1,ivar) = eu_i.*dCij.*ev_j;
            end              
        end
        
        function solvePDEs(obj)
            obj.adjointProb.setC(obj.homogenizedVariablesComputer.C);
            obj.adjointProb.computeVariables();
            obj.physicalProblem.setC(obj.homogenizedVariablesComputer.C);
            obj.physicalProblem.computeVariables();
        end
        
    end
    
    methods (Access = private)
        
        function createAdjointProblem(obj,fileName)
            fAdj = Preprocess.getBC_adjoint(fileName);
            obj.adjointProb = FEM.create(fileName);
            [dof,dofVal] = obj.adjointProb.dof.get_dof_conditions(fAdj,obj.adjointProb.dof.nunkn);
            obj.adjointProb.dof.neumann = dof;
            obj.adjointProb.dof.neumann_values = -dofVal;
            obj.adjointProb.preProcess;
        end
        
    end
end