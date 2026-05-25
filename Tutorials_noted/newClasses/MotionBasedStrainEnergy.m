classdef MotionBasedStrainEnergy < handle
    % Class to compute the motion based the strain energy and its gradient for a given
    % MP to then evaluate and constrain the motion based characterstic stiffness K=2*E/gamma^2
    % in MotionBasedStiffnessConstraint

    properties (Access = public)
        value0;
    end

    properties (Access = private)
        mesh
        filter
        material
        stateProblem % Elastic problem with prescribed displ (motion based)
        quadrature
        gradientFilter
        oldCost          
        oldGradient      
        xOld
    end

    methods (Access = public)
        function obj = MotionBasedStrainEnergy(cParams)
            obj.init(cParams);
            obj.createQuadrature();
        end

        function [J, dJ] = computeFunctionAndGradient(obj,x)
            xD = x.obtainDomainFunction();

            % Depending on the type of output of obtainDomainFunction
            if iscell(xD)
                xD_single = xD{1};
            else
                xD_single = xD;
            end

            % Filter the design variable
            x_filtered = obj.filter.compute(xD_single,2);
            xR = {x_filtered}; % convert it to cell for computeTensorAndGradient

            dx = xR{1} - obj.xOld;
            if norm(dx.fValues)/norm(xR{1}.fValues) > 0.05
                [C,dC] = obj.computeTensorFunctionAndGradient(xR);
                uS     = obj.computeStateVariable(C);
                J      = obj.computeFunctionValue(C,uS);
                dJ     = obj.computeGradient(dC{1},uS);
                dJ     = obj.gradientFilter.compute(dJ,2);
                dJval  = 0.5*dJ.fValues;
                if ~isempty(obj.value0)
                    dJval = dJval/obj.value0;
                end
                dJ.setFValues(dJval);
                dJ = {dJ};

                obj.oldCost     = J;
                obj.oldGradient = dJ;
                obj.xOld        = xR{1};
            else
                sp = ScalarProduct(obj.oldGradient{1}, dx, 'L2');
                J  = obj.oldCost + sp;
                dJ = obj.oldGradient;
            end
        end

        function title = getTitleToPlot(obj)
            title = 'Motion Based Strain Energy';
        end

    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh         = cParams.mesh;
            obj.filter       = cParams.filter;
            obj.material     = cParams.material;
            obj.stateProblem = cParams.stateProblem;
            obj.gradientFilter = cParams.gradientFilter;

            obj.xOld = 1000;

            % not needed anymore, it was to solve an error:
            if isfield(cParams,'value0')
                obj.value0 = cParams.value0;
            end
        end

        function createQuadrature(obj)
            quad = Quadrature.create(obj.mesh, 2);
            obj.quadrature = quad;
        end

        function xR = filterDesignVariable(obj, x)
            nDesVar = length(x);
            xR      = cell(nDesVar,1);
            for i = 1:nDesVar
                xR{i} = obj.filter.compute(x{i},2);
            end
        end

        function [C, dC] = computeTensorFunctionAndGradient(obj, xR)
            obj.material.setDesignVariable(xR);
            C  = obj.material.obtainTensor();
            dC = obj.material.obtainTensorDerivative();
        end
        
        function u = computeStateVariable(obj, C)
            obj.stateProblem.updateMaterial(C);
            obj.stateProblem.solve();
            u = obj.stateProblem.uFun;
        end
        
        function J = computeFunctionValue(obj, C, uS)
            stateStrain   = SymGrad(uS); 
            stress        = DDP(C, stateStrain); 
            dCompliance = DDP(stateStrain, stress); % corresponds to the amount of energy stored at each point of the domain
            
            % E = 1/2 * u * K * u
            J_abs = 0.5 * Integrator.compute(dCompliance, obj.mesh, obj.quadrature.order); % integrate using gauss quadrature (Koppen does it with matrices)
            
            if isempty(obj.value0)
                obj.value0 = J_abs;
            end
            J = J_abs / obj.value0;
        end

    end

    methods (Static, Access = private)
        function dJ = computeGradient(dC,uS)
            stateStrain   = SymGrad(uS);
            dStress       = DDP(dC,stateStrain);
            dJ            = DDP(stateStrain,dStress);
        end
    end



end


