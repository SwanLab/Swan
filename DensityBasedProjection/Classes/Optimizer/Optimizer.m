classdef Optimizer < handle

    properties (Access = private)
        mesh
        structure        
        projectorParameters
        filterParameters
        cost
        solverParameters
        iterations
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
                
                

                v        = sum(sum(obj.projectedField.D))/(obj.projectedField.volfracD*obj.mesh.elementNumberX*obj.mesh.elementNumberY)-1;
                dv_dxs   = 1/(obj.projectedField.volfracD*obj.mesh.elementNumberX*obj.mesh.elementNumberY)*ones(obj.mesh.elementNumberY,obj.mesh.elementNumberX);

                %   Derivate the proyected field respected the filteredField
                s.beta = obj.projectorParameters.beta;
                s.filteredField = obj.filteredField;

                s.eta =obj.projectorParameters.eta.E;
                E = ProjectedFieldDerivator(s);
                E.compute();
                derivatedProjectedField.E = E.derivatedProjectedField;

                s.eta =obj.projectorParameters.eta.I;
                I = ProjectedFieldDerivator(s);
                I.compute();
                derivatedProjectedField.I = I.derivatedProjectedField;

                s.eta =obj.projectorParameters.eta.D;
                D = ProjectedFieldDerivator(s);
                D.compute();
                derivatedProjectedField.D = D.derivatedProjectedField;

                
                %Derivate the filtered field by the field
                s.H = obj.filterParameters.H;
                s.Hs = obj.filterParameters.Hs;
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




                dvs_dxs(:)  = obj.filterParameters.H*(dv_dxs(:).*derivatedProjectedField.D(:)./obj.filterParameters.Hs);

                % Compute the solver
                sC = [];
                c = CostAndConstraintArturo(sC);
                [f0val,df0dx,df0dx2] = c.computeCostValueAndGradient(obj.mesh.elementNumberX,obj.mesh.elementNumberY);
                [fval,dfdx,dfdx2] =  c.computeConstraintValueAndGradient(obj.cost.E,obj.cost.I,obj.cost.D,v,penalizedDerivative.E,penalizedDerivative.I,penalizedDerivative.D,dvs_dxs,obj.cost.initial);
                xval        = obj.field(:);
                [xmma,~,zmma,~,~,~,~,~,~,obj.solverParameters.low,obj.solverParameters.upp] = mmasub2Arturo(obj.solverParameters.numberRestriction,obj.solverParameters.variableNumber,iter,xval,obj.solverParameters.minDensity,obj.solverParameters.maxDensity,obj.solverParameters.xold1,obj.solverParameters.xold2, ...
                    f0val,df0dx,df0dx2,fval,dfdx,dfdx2,obj.solverParameters.low,obj.solverParameters.upp,obj.solverParameters.a0,obj.solverParameters.mmaParameter.a,obj.solverParameters.mmaParameter.e,obj.solverParameters.mmaParameter.d);


                obj.solverParameters.xold2  = obj.solverParameters.xold1;
                obj.solverParameters.xold1  = xval;


                obj.field     = reshape(xmma,obj.mesh.elementNumberY,obj.mesh.elementNumberX);

                obj.cost.change = norm(xmma-xold,inf);

                obj.filteredField(:) = (obj.filterParameters.H*obj.field(:))./obj.filterParameters.Hs;

                %Project the new filtered field
             
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






                subplot(3,1,1);
                imagesc(-obj.projectedField.E); colormap(gray); caxis([-1 0]); axis off; axis equal; axis tight;
                subplot(3,1,2);
                imagesc(-obj.projectedField.I); colormap(gray); caxis([-1 0]); axis off; axis equal; axis tight;
                subplot(3,1,3);
                imagesc(-obj.projectedField.D); colormap(gray); caxis([-1 0]); axis off; axis equal; axis tight; pause(0.1)

                disp([' It.: ' sprintf('%4i',iter) ' Obj.: ' sprintf('%6.4f',zmma) ...
                    ' ui: '  sprintf('%12f', [obj.cost.E obj.cost.I obj.cost.D])...
                    ' V: '   sprintf('%6.3f',sum(obj.projectedField.I(:))/(obj.mesh.elementNumberX*obj.mesh.elementNumberY)) ...
                    ' ch.: ' sprintf('%6.3f', obj.cost.change)])

                if (obj.projectorParameters.beta < 32) && ((iterbeta >= 50) || (obj.cost.change <= 0.01))
                    obj.projectorParameters.beta     = 2*obj.projectorParameters.beta;
                    iterbeta = 0;
                    obj.cost.change   = 1;
                end

                if (itervol >= 20 )
                    itervol = 0;
                    obj.projectedField.volfracD = volumenFraction*sum(obj.projectedField.D(:))/sum(obj.projectedField.I(:));
                end
                if iter == obj.iterations               
                    finish = true;
                end

            end

        end 
    end 
end