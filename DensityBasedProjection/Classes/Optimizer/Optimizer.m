classdef Optimizer < handle

    properties (Access = private)
        mesh
        structure        
        projector
        filter
        cost
        solver
        iterations
        costChange
        field
        filteredField
        problemParameters
        
    end
    properties (Access = public)
        projectedField
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
                obj.mesh.elementNumberX = cParams.mesh.elementNumberX;
                obj.mesh.elementNumberY = cParams.mesh.elementNumberY;
                obj.mesh.neumanCondition = cParams.mesh.neumanCondition;
                obj.mesh.output = cParams.mesh.output;
                obj.mesh.freeDegress = cParams.mesh.freeDegress;
                obj.mesh.conectivityMatrixMat = cParams.mesh.conectivityMatrixMat;    

                obj.structure.elementType = cParams.structure.elementType;
                obj.structure.t=cParams.structure.t;
                obj.structure.penalization=cParams.structure.penalization;
                obj.structure.poissonCoefficient=cParams.structure.poissonCoefficient;
                obj.structure.elasticModuleMinimun=cParams.structure.elasticModuleMinimun;
                obj.structure.elasticModuleNeutral=cParams.structure.elasticModuleNeutral;
                obj.structure.elementalStiffnessMatrix=cParams.structure.elementalStiffnessMatrix;

                obj.projector.eta.E = cParams.projector.eta.E;
                obj.projector.eta.I = cParams.projector.eta.I;
                obj.projector.eta.D = cParams.projector.eta.D;
                obj.projector.beta = cParams.projector.beta;

                obj.filter.H = cParams.filter.H;
                obj.filter.Hs = cParams.filter.Hs;

                obj.cost = cParams.cost;
                obj.costChange = cParams.costChange;

                obj.solver.minDensity = cParams.solver.minDensity;
                obj.solver.maxDensity = cParams.solver.maxDensity;
                obj.solver.mmaParameter = cParams.solver.mmaParameter;
                obj.solver.numberRestriction = cParams.solver.numberRestriction;
                obj.solver.variableNumber = cParams.solver.variableNumber;
                obj.solver.xold1 = cParams.solver.xold1;
                obj.solver.xold2 = cParams.solver.xold2;
                obj.solver.low     = cParams.solver.low;
                obj.solver.upp     = cParams.solver.upp;
                obj.solver.a0      = cParams.solver.a0;

                obj.iterations = cParams.iterations;
                
                obj.field = cParams.field;
                obj.filteredField = cParams.filteredField;
                obj.projectedField= cParams.projectedField;

                obj.problemParameters.volfracD = cParams.problemParameters.volfracD;
        end 
        function optimize(obj)
            %Optimizer parameter
            iterbeta = 0;
            iter     = 0;
            itervol  = 0;
            finish = false;

           while (obj.costChange > 0.001) && (iter < 300) && (finish == false)
                %Actualización de variables al comienzo del bucle:
                iter     = iter     + 1;
                iterbeta = iterbeta + 1;
                itervol  = itervol  + 1;

                xold = obj.field(:);


                %Get cost and displacement
                s.mesh = obj.mesh;                                         
                s.structure= obj.structure;
                s.projectedField = obj.projectedField.E;   

                E = FEMcomputer(s);
                E.compute();
                obj.cost.E = E.cost;
                displacement.E = E.displacement;

                s.projectedField = obj.projectedField.I;                      
                I = FEMcomputer(s);
                I.compute();
                obj.cost.I = I.cost;
                displacement.I = I.displacement;

                s.projectedField = obj.projectedField.D;                      
                D = FEMcomputer(s);
                D.compute();
                obj.cost.D = D.cost;
                displacement.D = D.displacement;
                           
              %Get the cost derivated respective the proyectedField 

                s.elementNumberX = obj.mesh.elementNumberX;
                s.elementNumberY = obj.mesh.elementNumberY;
                s.elementalStiffnessMatrix = obj.structure.elementalStiffnessMatrix;
                s.conectivityMatrixMat = obj.mesh.conectivityMatrixMat; 

                s.displacement =displacement.E;
                E = ProjectedCostDerivator(s);
                E.compute();
                derivedCost.E = E.derivedCost;

                s.displacement =displacement.I;
                I = ProjectedCostDerivator(s);
                I.compute();
                derivedCost.I = I.derivedCost;

                s.displacement =displacement.D;
                D = ProjectedCostDerivator(s);
                D.compute();
                derivedCost.D = D.derivedCost;

                % Penalize the derivative cost respect the proyectedField 
               
                s.elasticModuleMinimun = obj.structure.elasticModuleMinimun;
                s.elasticModuleNeutral = obj.structure.elasticModuleNeutral;
                s.penalization = obj.structure.penalization;

                s.nonPenalizedVariable = derivedCost.E;
                s.projectedField = obj.projectedField.E ;
                E = DerivativePenalizer(s);
                E.penalize();
                penalizedDerivative.E = E.penalizedDerivative;

                s.nonPenalizedVariable = derivedCost.I;
                s.projectedField = obj.projectedField.I ;
                I = DerivativePenalizer(s);
                I.penalize();
                penalizedDerivative.I = I.penalizedDerivative;

                s.nonPenalizedVariable = derivedCost.D;
                s.projectedField = obj.projectedField.D ;
                D = DerivativePenalizer(s);
                D.penalize();
                penalizedDerivative.D = D.penalizedDerivative;
                
                

                v        = sum(sum(obj.projectedField.D))/(obj.problemParameters.volfracD*obj.mesh.elementNumberX*obj.mesh.elementNumberY)-1;
                dv_dxs   = 1/(obj.problemParameters.volfracD*obj.mesh.elementNumberX*obj.mesh.elementNumberY)*ones(obj.mesh.elementNumberY,obj.mesh.elementNumberX);

                %   Derivate the proyected field respected the filteredField
                s.beta = obj.projector.beta;
                s.filteredField = obj.filteredField;

                s.eta =obj.projector.eta.E;
                E = ProjectedFieldDerivator(s);
                E.compute();
                derivatedProjectedField.E = E.derivatedProjectedField;

                s.eta =obj.projector.eta.I;
                I = ProjectedFieldDerivator(s);
                I.compute();
                derivatedProjectedField.I = I.derivatedProjectedField;

                s.eta =obj.projector.eta.D;
                D = ProjectedFieldDerivator(s);
                D.compute();
                derivatedProjectedField.D = D.derivatedProjectedField;

                
                %Derivate the filtered field by the field
                s.H = obj.filter.H;
                s.Hs = obj.filter.Hs;
                B = FilteredFieldDerivator(s);
                B.compute();
                derivedFilteredField = B.derivedFilteredField;
              
                % Calculate the cost derivated by the field
                s.derivedFilteredField = derivedFilteredField;

                s.derivedProyectedCost = penalizedDerivative.E;
                s.derivatedProjectedField = derivatedProjectedField.E;                
                E = CostDerivator(s);
                E.compute();
                penalizedDerivative.E = E.derivedCost;

                s.derivedProyectedCost = penalizedDerivative.I;
                s.derivatedProjectedField = derivatedProjectedField.I;                
                I = CostDerivator(s);
                I.compute();
                penalizedDerivative.I = I.derivedCost;

                s.derivedProyectedCost = penalizedDerivative.D;
                s.derivatedProjectedField = derivatedProjectedField.D;                
                D = CostDerivator(s);
                D.compute();
                penalizedDerivative.D = D.derivedCost;




                dvs_dxs(:)  = obj.filter.H*(dv_dxs(:).*derivatedProjectedField.D(:)./obj.filter.Hs);

                % Compute the solver
                sC = [];
                c = CostAndConstraintArturo(sC);
                [f0val,df0dx,df0dx2] = c.computeCostValueAndGradient(obj.mesh.elementNumberX,obj.mesh.elementNumberY);
                [fval,dfdx,dfdx2] =  c.computeConstraintValueAndGradient(obj.cost.E,obj.cost.I,obj.cost.D,v,penalizedDerivative.E,penalizedDerivative.I,penalizedDerivative.D,dvs_dxs,obj.cost.initial);
                xval        = obj.field(:);
                [xmma,~,zmma,~,~,~,~,~,~,obj.solver.low,obj.solver.upp] = mmasub2Arturo(obj.solver.numberRestriction,obj.solver.variableNumber,iter,xval,obj.solver.minDensity,obj.solver.maxDensity,obj.solver.xold1,obj.solver.xold2, ...
                    f0val,df0dx,df0dx2,fval,dfdx,dfdx2,obj.solver.low,obj.solver.upp,obj.solver.a0,obj.solver.mmaParameter.a,obj.solver.mmaParameter.e,obj.solver.mmaParameter.d);


                obj.solver.xold2  = obj.solver.xold1;
                obj.solver.xold1  = xval;


                obj.field     = reshape(xmma,obj.mesh.elementNumberY,obj.mesh.elementNumberX);

                obj.costChange = norm(xmma-xold,inf);

                obj.filteredField(:) = (obj.filter.H*obj.field(:))./obj.filter.Hs;

                %Project the new filtered field
             
                s.beta = obj.projector.beta;
                s.eta = obj.projector.eta.E;
                s.filteredField =obj.filteredField;
                E = FieldProjector(s);
                E.compute();
                obj.projectedField.E = E.projectedField;
    
                s.beta =  obj.projector.beta;
                s.eta = obj.projector.eta.I;
                I = FieldProjector(s);
                I.compute();
                obj.projectedField.I = I.projectedField;
    
                s.beta =  obj.projector.beta;
                s.eta = obj.projector.eta.D;
                D = FieldProjector(s);
                D.compute();
                obj.projectedField.D = D.projectedField;






                subplot(3,1,1);
                imagesc(-obj.projectedField.E); colormap(gray); caxis([-1 0]); axis off; axis equal; axis tight;
                subplot(3,1,2);
                imagesc(-obj.projectedField.I); colormap(gray); caxis([-1 0]); axis off; axis equal; axis tight;
                subplot(3,1,3);
                imagesc(-obj.projectedField.D); colormap(gray); caxis([-1 0]); axis off; axis equal; axis tight; pause(0.1)

                disp([' It.: ' sprintf('%4i',iter) ' Obj.: ' sprintf('%6.4f',zmma) ...
                    ' ui: '  sprintf('%12f', [obj.cost.E obj.cost.I obj.cost.D])...
                    ' V: '   sprintf('%6.3f',sum(obj.projectedField.I(:))/(obj.mesh.elementNumberX*obj.mesh.elementNumberY)) ...
                    ' ch.: ' sprintf('%6.3f', obj.costChange)])

                if (obj.projector.beta < 32) && ((iterbeta >= 50) || (obj.costChange <= 0.01))
                    obj.projector.beta     = 2*obj.projector.beta;
                    iterbeta = 0;
                    obj.costChange   = 1;
                end

                if (itervol >= 20 )
                    itervol = 0;
                    obj.problemParameters.volfracD = volumenFraction*sum(obj.projectedField.D(:))/sum(obj.projectedField.I(:));
                end
                if iter == obj.iterations               
                    finish = true;
                end

            end

        end 
    end 
end