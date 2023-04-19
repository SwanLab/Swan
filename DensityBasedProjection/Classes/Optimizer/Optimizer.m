classdef Optimizer < handle
    properties (Access = public)
        E
        I
        D
    end

    properties (Access = private)
        mesh
        structure
        projectorParameters
        filterParameters
        displacement
        solverParameters
        iterations
        problemParameters
        newField
    end

    methods (Access = public)
        function obj = Optimizer(cParams)
            obj.inputData(cParams);
        end

        function compute(obj)
            obj.optimize();
        end
    end
    methods (Access = private)
        function inputData(obj,cParams)

            obj.mesh = cParams.mesh;
            obj.structure = cParams.structure;
            obj.structure.elementType = 'Square';
            obj.projectorParameters = cParams.projectorParameters;
            obj.filterParameters= cParams.filterParameters;
            obj.solverParameters = cParams.solverParameters;
            obj.iterations = cParams.iterations;
            obj.E = cParams.E;
            obj.I = cParams.I;
            obj.D = cParams.D;
        end
        function optimize(obj)
            %Optimizer parameter
            iterbeta = 0;
            iter     = 0;
            itervol  = 0;
            finish = false;

            while (obj.solverParameters.costChange > 0.001) && (iter < 300) && (finish == false)
                %Actualización de variables al comienzo del bucle:
                iter     = iter     + 1;
                iterbeta = iterbeta + 1;
                itervol  = itervol  + 1;

                oldE = obj.E;
                oldI = obj.I;
                oldD = obj.D;

                % Compute cost and displacement with the projected field
                obj.computeCost();
                % Derivate the projected field respected the filteredField
                obj.deriveProjectedFieldRespectedFilteredField();
                % Calculate the cost derivated by the field
                obj.deriveCostRespectedField();
                % Compute volumen, derivate volumen, filter derivated volumen
                obj.computeVolumenValues();
                % Compute the solver
                obj.computeSolver(iter,oldE.designField.field);
                %Set the new field values
                obj.setNewFields();
                %Filter the new field
                obj.filterNewField();
                %Project the new filtered field
                obj.projectNewFilteredField();
                %Set the design Field in the design objects
                obj.updateDataInDesignVariables();

                %Plot results
                obj.plotResults(iter);
                %New optimizer parameters
                [iterbeta,itervol,finish,iter] = Optimizer.reconfigurateOptimizeParameters(obj,iterbeta,itervol,finish,iter);


            end

        end
        function updateDataInDesignVariables(obj)

            obj.E.designCost.designFields = obj.E.designField;
            obj.I.designCost.designFields = obj.I.designField;
            obj.D.designCost.designFields = obj.D.designField;

            obj.E.designVolumen.designFields = obj.E.designField;
            obj.I.designVolumen.designFields = obj.I.designField;
            obj.D.designVolumen.designFields = obj.D.designField;

        end

        function computeCost(obj)
            %Get cost and displacement
            obj.E.designCost.computeCost();
            obj.I.designCost.computeCost();
            obj.D.designCost.computeCost();
        end

        function deriveCostRespectedField(obj)
            obj.E.designCost.deriveCost();
            obj.I.designCost.deriveCost();
            obj.D.designCost.deriveCost();
        end
        function computeVolumenValues(obj)
            obj.E.designVolumen.computeVolumen();
            obj.I.designVolumen.computeVolumen();
            obj.D.designVolumen.computeVolumen();

            obj.E.designVolumen.deriveVolumenAndFilter();
            obj.I.designVolumen.deriveVolumenAndFilter();
            obj.D.designVolumen.deriveVolumenAndFilter();

        end
        function computeSolver(obj,iter,xold)
            sC = [];
            c = SolverMMA(sC);
            [f0val,df0dx,df0dx2] = c.computeCostValueAndGradient(obj.mesh.elementNumberX,obj.mesh.elementNumberY);
            [fval,dfdx,dfdx2] =  c.computeConstraintValueAndGradient(obj.E.designCost.cost,obj.I.designCost.cost,obj.D.designCost.cost,obj.D.designVolumen.volumen,obj.E.designCost.derivedCost,obj.I.designCost.derivedCost,obj.D.designCost.derivedCost,obj.D.designVolumen.derivedVolumen,obj.solverParameters.initialCost);

            xval        = obj.E.designField.field(:);
            [obj.solverParameters.xmma,~,obj.solverParameters.zmma,~,~,~,~,~,~,obj.solverParameters.low,obj.solverParameters.upp] = mmasub2Arturo(obj.solverParameters.numberRestriction,obj.solverParameters.variableNumber,iter,xval,obj.solverParameters.minDensity,obj.solverParameters.maxDensity,obj.solverParameters.xold1,obj.solverParameters.xold2, ...
                f0val,df0dx,df0dx2,fval,dfdx,dfdx2,obj.solverParameters.low,obj.solverParameters.upp,obj.solverParameters.a0,obj.solverParameters.mmaParameter.a,obj.solverParameters.mmaParameter.e,obj.solverParameters.mmaParameter.d);


            obj.solverParameters.xold2  = obj.solverParameters.xold1;
            obj.solverParameters.xold1  = xval;
            obj.solverParameters.costChange  = norm(obj.solverParameters.xmma-xold(:),inf);

            obj.newField     = reshape(obj.solverParameters.xmma,obj.mesh.elementNumberY,obj.mesh.elementNumberX);
        end
        function setNewFields(obj)
            obj.E.designField.field = obj.newField;
            obj.I.designField.field = obj.newField;
            obj.D.designField.field = obj.newField;
        end
        function filterNewField(obj)
            obj.E.designField.filter();
            obj.I.designField.filter();
            obj.D.designField.filter();
        end
        function projectNewFilteredField(obj)
            obj.E.designField.project();
            obj.I.designField.project();
            obj.D.designField.project();

        end
        function deriveProjectedFieldRespectedFilteredField(obj)
            obj.E.designField.deriveProjectedField();
            obj.I.designField.deriveProjectedField();
            obj.D.designField.deriveProjectedField();
        end


        function plotResults(obj,iter)
            %             subplot(3,1,1);
            %             imagesc(-obj.projectedField.E); colormap(gray); caxis([-1 0]); axis off; axis equal; axis tight;
            %             subplot(3,1,2);
            %             imagesc(-obj.projectedField.I); colormap(gray); caxis([-1 0]); axis off; axis equal; axis tight;
            %             subplot(3,1,3);
            %             imagesc(-obj.projectedField.D); colormap(gray); caxis([-1 0]); axis off; axis equal; axis tight; pause(0.1)
            disp([' It.: ' sprintf('%4i',iter) ' Obj.: ' sprintf('%6.4f',obj.solverParameters.zmma) ...
                ' ui: '  sprintf('%12f', [obj.E.designCost.cost obj.I.designCost.cost obj.D.designCost.cost])...
                ' V: '   sprintf('%6.3f',sum(obj.I.designField.projectedField(:))/(obj.mesh.elementNumberX*obj.mesh.elementNumberY)) ...
                ' ch.: ' sprintf('%6.3f', obj.solverParameters.costChange)])
        end
    end
    methods (Static)
        function [iterbeta,itervol,finish,iter] = reconfigurateOptimizeParameters(obj,iterbeta,itervol,finish,iter)
            if (obj.projectorParameters.beta < 32) && ((iterbeta >= 50) || (obj.solverParameters.costChange <= 0.01))
                obj.projectorParameters.beta     = 2*obj.projectorParameters.beta;
                iterbeta = 0;
                obj.solverParameters.costChange   = 1;
            end

            if (itervol >= 20 )
                itervol = 0;
                E.designVolumen.computeVolumenFraction(D,I);
                I.designVolumen.computeVolumenFraction(D,I);
                D.designVolumen.computeVolumenFraction(D,I);
            end
            if iter == obj.iterations
                finish = true;
            end
        end
    end
end