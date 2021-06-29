classdef ShFunc_Chomog < ShapeFunctional
    
    properties (Access = public)
        Chomog
        tstress
        tstrain
        Chomog_Derivatives
        physicalProblem
        rho
        matProps
    end
    
    properties (Access = protected)
        regDesignVariable
    end    
    
    methods (Access = public)
        
        
        function f = getPhysicalProblems(obj)
            f{1} = obj.physicalProblem;            
        end
        
        function f = getRegularizedDensity(obj)
            f = obj.regDesignVariable;
        end
        
        function fP = addPrintableVariables(obj)
            fP = obj.addMicroVariables();
            fP{end+1}.value = obj.gradient;
            fH = obj.homogenizedVariablesComputer.addPrintableVariables(obj.designVariable);
            for i = 1:numel(fH)
                fP{end+1} = fH{i};
            end
        end
        
        function fP = addMicroVariables(obj)           
            p = obj.getPhysicalProblems();
            p = p{1};
            nMicro = numel(p.variables2print);
            fP = cell(nMicro,1);
            for i = 1:nMicro
                pV = p;
                pV.variables = p.variables2print{i};
                cParams.physicalProblem = pV;
                g = PdeVariableToPrintGetter(cParams);
                v = g.compute();
                fP{i}.value = v;            
            end
        end
        
        function fP = createPrintVariables(obj)           
            types = {'Elasticity','Elasticity','Elasticity','ScalarNodal'};      
            names = {'Primal1','Primal2','Primal3','Gradient'};
            
            fP = obj.obtainPrintVariables(types,names); 
            fP = obj.addHomogPrintVariablesNames(fP);            
        end   
        
        function v = getVariablesToPlot(obj)
            v{1} = obj.value*obj.value0;
        end

        function t = getTitlesToPlot(obj)
            t{1} = 'C - C not scaled';
        end    
        
        function q = getQuad(obj)
            q = obj.physicalProblem.element.quadrature;
        end        
        
    end
    
    methods (Access = protected)
        
        
        function compute_Chomog_Derivatives(obj)
             dC = obj.homogenizedVariablesComputer.dC;
            
            nstre = obj.physicalProblem.element.getNstre();
            nelem = obj.physicalProblem.mesh.nelem;
            ngaus = obj.physicalProblem.element.quadrature.ngaus;
            
            
            obj.Chomog_Derivatives = zeros(nstre,nstre,ngaus,nelem);
            for istreChomog = 1:nstre
                for jstreChomog = 1:nstre
                    for igaus=1:ngaus
                        for istre=1:nstre
                            for jstre = 1:nstre
                                obj.Chomog_Derivatives(istreChomog,jstreChomog,igaus,:) = ...
                                    squeeze(obj.Chomog_Derivatives(istreChomog,jstreChomog,igaus,:)) + ...
                                    (squeeze(obj.tstrain(istreChomog,igaus,istre,:))...
                                    .*squeeze(dC(istre,jstre,:,:,igaus))...
                                    .*squeeze(obj.tstrain(jstreChomog,igaus,jstre,:)));
                            end
                        end
                    end
                end
            end
            
        end
        
        function initChomog(obj,cParams)
            cParams.filterParams.quadratureOrder = 'LINEAR';
            obj.init(cParams);
            obj.physicalProblem = FEM.create(cParams.femSettings.fileName);            
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
        
        function computePhysicalData(obj)
            obj.updateHomogenizedMaterialProperties();
            obj.solveState();
        end
        
        function solveState(obj)
            obj.physicalProblem.setC(obj.homogenizedVariablesComputer.C)
            obj.physicalProblem.computeChomog();
            obj.Chomog  = obj.physicalProblem.variables.Chomog;
            obj.tstrain = obj.physicalProblem.variables.tstrain;
            obj.tstress = obj.physicalProblem.variables.tstress;            
        end
        
        function updateHomogenizedMaterialProperties(obj)
            rhoV{1} = obj.filter.getP0fromP1(obj.designVariable.value);
            obj.regDesignVariable = rhoV{1};           
            obj.homogenizedVariablesComputer.computeCtensor(rhoV);
            obj.homogenizedVariablesComputer.computeDensity(rhoV);            
        end
        
        
   
    end
    
    methods (Access = protected, Static)
        function r = projection_Chomog(inv_matCh,alpha,beta)
            weights = alpha*beta';
            r = sum(sum(weights.*inv_matCh));
        end
    end
end
