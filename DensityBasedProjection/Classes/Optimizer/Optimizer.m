classdef Optimizer < handle
    properties (Access = public)
        projectedField
    end

    properties (Access = private)
        mesh
        structure
        projectorParameters
        filterParameters
        cost
        derivedCost
        displacement
        solverParameters
        iterations
        field
        filteredField
        problemParameters
        derivedProjectedField
        derivedFilteredField
        volumen

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
            obj.projectorParameters = cParams.projector;
            obj.filterParameters= cParams.filterParameters;
            obj.cost = cParams.cost;
            obj.solverParameters = cParams.solverParameters;
            obj.iterations = cParams.iterations;
            obj.field = cParams.field;
            obj.filteredField = cParams.filteredField;
            obj.projectedField= cParams.projectedField;
            obj.volumen = cParams.volumen;
        end
        function optimize(obj)
            %Optimizer parameter
            iterbeta = 0;
            iter     = 0;
            itervol  = 0;
            finish = false;

            while (obj.cost.change > 0.001) && (iter < 300) && (finish == false)
                %Actualización de variables al comienzo del bucle:
                iter     = iter     + 1;
                iterbeta = iterbeta + 1;
                itervol  = itervol  + 1;

                xold = obj.field(:);

                % Compute cost and displacement with the projected field
                obj.computeFEM();
                % Derivate the projected field respected the filteredField
                obj.deriveProjectedFieldRespectedFilteredField();
                % Calculate the cost derivated by the field
                obj.deriveCostRespectedField();
                % Compute volumen, derivate volumen, filter derivated volumen
                obj.computeVolumenValues();

                % Compute the solver
                obj.computeSolver(iter,xold);

                %Filter the new field
                obj.filterNewField();
                %Project the new filtered field
                obj.projectNewFilteredField();
                %Plot results
                obj.plotResults(iter);
                %New optimizer parameters
                [iterbeta,itervol,finish,iter] = Optimizer.reconfigurateOptimizeParameters(obj,iterbeta,itervol,finish,iter);

            end

        end
    
        function computeFEM(obj)
            %Get cost and displacement
            s.mesh = obj.mesh;
            s.structure= obj.structure;
            s.projectedField = obj.projectedField.E;

            E = FEMcomputer(s);
            E.compute();
            obj.cost.E = E.cost;
            obj.displacement.E = E.displacement;

            s.projectedField = obj.projectedField.I;
            I = FEMcomputer(s);
            I.compute();
            obj.cost.I = I.cost;
            obj.displacement.I = I.displacement;

            s.projectedField = obj.projectedField.D;
            D = FEMcomputer(s);
            D.compute();
            obj.cost.D = D.cost;
            obj.displacement.D = D.displacement;

        end

        function deriveCostRespectedField(obj)
            s.mesh = obj.mesh;
            s.structure =obj.structure;
            s.filterParameters =obj.filterParameters ;
            s.filteredField =obj.filteredField ;
            s.projectorParameters.beta =obj.projectorParameters.beta;
            s.cost =obj.cost ;

            
            s.displacement =obj.displacement.E;
            s.projectedField =obj.projectedField.E;
            s.derivedProjectedField = obj.derivedProjectedField.E;
            s.projectorParameters.eta =obj.projectorParameters.eta.E ;

            E = CostFieldDerivator(s);
            E.compute();
            obj.derivedCost.E = E.derivedCost;

            s.displacement =obj.displacement.I;
            s.projectedField =obj.projectedField.I;
            s.derivedProjectedField = obj.derivedProjectedField.I;
            s.projectorParameters.eta =obj.projectorParameters.eta.I ;
            I = CostFieldDerivator(s);
            I.compute();
            obj.derivedCost.I = I.derivedCost;

            s.displacement =obj.displacement.D;
            s.projectedField =obj.projectedField.D;
            s.projectorParameters.eta =obj.projectorParameters.eta.D;
            s.derivedProjectedField = obj.derivedProjectedField.D;
            D = CostFieldDerivator(s);
            D.compute();
            obj.derivedCost.D = D.derivedCost;

        end
        function computeSolver(obj,iter,xold)
            sC = [];
            c = CostAndConstraintArturo(sC);
            [f0val,df0dx,df0dx2] = c.computeCostValueAndGradient(obj.mesh.elementNumberX,obj.mesh.elementNumberY);
            [fval,dfdx,dfdx2] =  c.computeConstraintValueAndGradient(obj.cost.E,obj.cost.I,obj.cost.D,obj.volumen.value,obj.derivedCost.E,obj.derivedCost.I,obj.derivedCost.D,obj.volumen.derivated,obj.cost.initial);
            xval        = obj.field(:);
            [obj.solverParameters.xmma,~,obj.solverParameters.zmma,~,~,~,~,~,~,obj.solverParameters.low,obj.solverParameters.upp] = mmasub2Arturo(obj.solverParameters.numberRestriction,obj.solverParameters.variableNumber,iter,xval,obj.solverParameters.minDensity,obj.solverParameters.maxDensity,obj.solverParameters.xold1,obj.solverParameters.xold2, ...
                f0val,df0dx,df0dx2,fval,dfdx,dfdx2,obj.solverParameters.low,obj.solverParameters.upp,obj.solverParameters.a0,obj.solverParameters.mmaParameter.a,obj.solverParameters.mmaParameter.e,obj.solverParameters.mmaParameter.d);


            obj.solverParameters.xold2  = obj.solverParameters.xold1;
            obj.solverParameters.xold1  = xval;
            obj.field     = reshape(obj.solverParameters.xmma,obj.mesh.elementNumberY,obj.mesh.elementNumberX);
            obj.cost.change = norm(obj.solverParameters.xmma-xold,inf);

        end
        function computeVolumenValues(obj)
            s.mesh =obj.mesh;
            s.filterParameters = obj.filterParameters;
            s.projectedField = obj.projectedField;
            s.derivedProjectedField = obj.derivedProjectedField;
            s.volumen = obj.volumen;
            B = VolumenComputer(s);
            B.compute;
            obj.volumen = B.volumen;
        end
        function deriveProjectedFieldRespectedFilteredField(obj)
            s.beta = obj.projectorParameters.beta;
            s.filteredField = obj.filteredField;

            s.eta =obj.projectorParameters.eta.E;
            E = ProjectedFieldFilteredFieldDerivator(s);
            E.compute();
            obj.derivedProjectedField.E = E.derivatedProjectedField;

            s.eta =obj.projectorParameters.eta.I;
            I = ProjectedFieldFilteredFieldDerivator(s);
            I.compute();
            obj.derivedProjectedField.I = I.derivatedProjectedField;

            s.eta =obj.projectorParameters.eta.D;
            D = ProjectedFieldFilteredFieldDerivator(s);
            D.compute();
            obj.derivedProjectedField.D = D.derivatedProjectedField;
        end
        function projectNewFilteredField(obj)
            s.beta = obj.projectorParameters.beta;
            s.eta = obj.projectorParameters.eta.E;
            s.filteredField =obj.filteredField;
            E = FieldProjector(s);
            E.compute();
            obj.projectedField.E = E.projectedField;

            s.beta =  obj.projectorParameters.beta;
            s.eta = obj.projectorParameters.eta.I;
            I = FieldProjector(s);
            I.compute();
            obj.projectedField.I = I.projectedField;

            s.beta =  obj.projectorParameters.beta;
            s.eta = obj.projectorParameters.eta.D;
            D = FieldProjector(s);
            D.compute();
            obj.projectedField.D = D.projectedField;

        end
        function filterNewField(obj)
            s.filterParameters = obj.filterParameters;
            s.field = obj.field;
            B = FilterComputer(s);
            B.compute();
            obj.filteredField = B.filteredField;
        end

        function plotResults(obj,iter)
%             subplot(3,1,1);
%             imagesc(-obj.projectedField.E); colormap(gray); caxis([-1 0]); axis off; axis equal; axis tight;
%             subplot(3,1,2);
%             imagesc(-obj.projectedField.I); colormap(gray); caxis([-1 0]); axis off; axis equal; axis tight;
%             subplot(3,1,3);
%             imagesc(-obj.projectedField.D); colormap(gray); caxis([-1 0]); axis off; axis equal; axis tight; pause(0.1)
            disp([' It.: ' sprintf('%4i',iter) ' Obj.: ' sprintf('%6.4f',obj.solverParameters.zmma) ...
                ' ui: '  sprintf('%12f', [obj.cost.E obj.cost.I obj.cost.D])...
                ' V: '   sprintf('%6.3f',sum(obj.projectedField.I(:))/(obj.mesh.elementNumberX*obj.mesh.elementNumberY)) ...
                ' ch.: ' sprintf('%6.3f', obj.cost.change)])
        end
    end
    methods (Static)
        function [iterbeta,itervol,finish,iter] = reconfigurateOptimizeParameters(obj,iterbeta,itervol,finish,iter)
            if (obj.projectorParameters.beta < 32) && ((iterbeta >= 50) || (obj.cost.change <= 0.01))
                obj.projectorParameters.beta     = 2*obj.projectorParameters.beta;
                iterbeta = 0;
                obj.cost.change   = 1;
            end

            if (itervol >= 20 )
                itervol = 0;
                obj.volumen.volfracD = volumenFraction*sum(obj.projectedField.D(:))/sum(obj.projectedField.I(:));
            end
            if iter == obj.iterations
                finish = true;
            end
        end
    end
end