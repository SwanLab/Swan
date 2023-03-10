classdef ComplianceRobustComputer < handle
    properties (Access = public)
        projectedField
        iterations
    end
    properties (Access = private)
    end
    methods (Access = public)
        function obj = ComplianceRobustComputer(cParams)
            obj.inputData(cParams);
        end
        function compute(obj)
            obj.runCode();
        end
    end
    methods (Access = private)
        function inputData(obj,cParams)
            obj.iterations = cParams.iterations;
        end 
        
        function runCode(obj)

            %% Parameters

            %Malla (elementos cuadrados) y dimensiones de la estructura:
            %Malla estandar 400x80
            elementNumberX    = 160;
            elementNumberY    = 80;
            t       = 1;

            %Parámetros del material:
            elasticModuleNeutral      = 1;
            elasticModuleMinimun    = 1e-6;
            poissonCoefficient      = 0.3;
            penalization   = 3;

            %Entradas del problema:
            volumenFraction =   0.5;
            minimunInfluenceRadios    =   1.5;
            neumanCondition       = -1e-3;
            output    = 2;

            %Inicializaciones: %Es una clase de desplazamientos; 2 dof por cada elemento
            displacementE = zeros(2*(elementNumberY+1)*(elementNumberX+1),1);
            displacementI = zeros(2*(elementNumberY+1)*(elementNumberX+1),1);
            displacementD = zeros(2*(elementNumberY+1)*(elementNumberX+1),1);
            field        = volumenFraction*ones(elementNumberY,elementNumberX);
            filteredField = field;

            %Define Field and Desity
            field([],[]) = 1;
            minDensity    = zeros(elementNumberY,elementNumberX);
            maxDensity    =  ones(elementNumberY,elementNumberX);
            minDensity = minDensity(:);
            maxDensity = maxDensity(:);

            %MMA parameters 
            costChange  = 1;
            numberRestriction       = 4;
            variableNumber       = elementNumberX*elementNumberY;
            low     = minDensity;
            upp     = maxDensity;
            a0      = 1;
            mmaParameter.a       = [1 1 1 0]';
            mmaParameter.e       = 1000*ones(numberRestriction,1);
            mmaParameter.d       = 0*ones(numberRestriction,1);
            xold1   = field(:);
            xold2   = field(:);

            %Elemental stifness Matrix
            s.elementType = 'square';
            s.t = t;
            s.poissonCoefficient = poissonCoefficient;
            B = ElementalStiffnessMatricesComputer(s);
            B.compute();
            elementalStiffnessMatrix = B.elementalStiffnessMatrix;

            %Mesh parameters
            s.elementNumberX =  elementNumberX;
            s.elementNumberY =  elementNumberY;
            B = GeometryComputer(s);
            B.compute();
            allDegress   = B.degress.all;
            freeDegress  = B.degress.free;
            fixedDegress = B.degress.fixed;
            conectivityMatrixMat = B.conectivityMatrixMat;


            %Filter and Proyecition parameters:
            s.elementNumberX =  elementNumberX;
            s.elementNumberY =  elementNumberY;
            s.minimunInfluenceRadios = minimunInfluenceRadios;
            B = weightFilterComputer(s);
            B.compute();
            H = B.H;
            Hs = B.Hs;

            eta      = 0.25;
            etaE     = 1-eta;
            etaI     = 0.5;
            etaD     = eta;
            beta     = 1;
            %Optimizer parameter
            iterbeta = 0;
            iter     = 0;
            itervol  = 0;
            finish = false;

            % Project the initial field
            s.beta = beta;
            s.eta =etaE;
            s.filteredField =filteredField;
            E = FieldProjector(s);
            E.compute();
            projectedFieldE = E.projectedField;

            s.beta = beta;
            s.eta =etaI;
            s.filteredField =filteredField;
            I = FieldProjector(s);
            I.compute();
            projectedFieldI = I.projectedField;

            s.beta = beta;
            s.eta =etaD;
            s.filteredField =filteredField;
            D = FieldProjector(s);
            D.compute();
            projectedFieldD = D.projectedField;

            volfracD = volumenFraction*sum(projectedFieldD(:))/sum(projectedFieldI(:));

            %Get intial cost
            s.mesh.elementNumberX = elementNumberX;
            s.mesh.elementNumberY = elementNumberY;
            s.mesh.neumanCondition = neumanCondition;
            s.mesh.output = output;
            s.mesh.freeDegress = freeDegress;
            s.mesh.conectivityMatrixMat = conectivityMatrixMat;
            
            s.structure.elementalStiffnessMatrix = elementalStiffnessMatrix;
            s.structure.t=t;
            s.structure.penalization=penalization;
            s.structure.poissonCoefficient=poissonCoefficient;
            s.structure.elasticModuleMinimun=elasticModuleMinimun;
            s.structure.elasticModuleNeutral=elasticModuleNeutral;
            s.projectedField = projectedFieldE;
            
            B = FEMcomputer(s);
            B.compute();
            cost = B.cost;
            

            %Optimización
            while (costChange > 0.001) && (iter < 300) && (finish == false)

                %Actualización de variables al comienzo del bucle:
                iter     = iter     + 1;
                iterbeta = iterbeta + 1;
                itervol  = itervol  + 1;

                xold = field(:);


             %Get cost and displacement
                s.mesh.elementNumberX = elementNumberX;
                s.mesh.elementNumberY = elementNumberY;
                s.mesh.neumanCondition = neumanCondition;
                s.mesh.output = output;
                s.mesh.freeDegress = freeDegress;
                s.mesh.conectivityMatrixMat = conectivityMatrixMat;                            
                s.structure.elementType = 'square';
                s.structure.t=t;
                s.structure.penalization=penalization;
                s.structure.poissonCoefficient=poissonCoefficient;
                s.structure.elasticModuleMinimun=elasticModuleMinimun;
                s.structure.elasticModuleNeutral=elasticModuleNeutral;
                
                s.projectedField = projectedFieldE;                      
                E = FEMcomputer(s);
                E.compute();
                costE = E.cost;
                displacementE = E.displacement;

                s.projectedField = projectedFieldI;                      
                I = FEMcomputer(s);
                I.compute();
                costI = I.cost;
                displacementI = I.displacement;

                s.projectedField = projectedFieldD;                      
                D = FEMcomputer(s);
                D.compute();
                costD = D.cost;
                displacementD = D.displacement;
                           
              %Get the cost derivated respective the proyectedField 

                s.elementNumberX = elementNumberX;
                s.elementNumberY = elementNumberY;
                s.elementalStiffnessMatrix = elementalStiffnessMatrix;
                s.conectivityMatrixMat = conectivityMatrixMat; 

                s.displacement =displacementE;
                E = ProjectedCostDerivator(s);
                E.compute();
                cE1a = E.derivedCost;

                s.displacement =displacementI;
                I = ProjectedCostDerivator(s);
                I.compute();
                cI1a = I.derivedCost;

                s.displacement =displacementD;
                D = ProjectedCostDerivator(s);
                D.compute();
                cD1a = D.derivedCost;

              % Penalize the derivative cost respect the proyectedField 
               
                s.elasticModuleMinimun = elasticModuleMinimun;
                s.elasticModuleNeutral = elasticModuleNeutral;
                s.penalization = penalization;

                s.nonPenalizedVariable = cE1a;
                s.projectedField = projectedFieldE ;
                E = DerivativePenalizer(s);
                E.penalize();
                dcsE1 = E.penalizedDerivative;

                s.nonPenalizedVariable = cI1a;
                s.projectedField = projectedFieldI ;
                I = DerivativePenalizer(s);
                I.penalize();
                dcsI1 = I.penalizedDerivative;

                s.nonPenalizedVariable = cD1a;
                s.projectedField = projectedFieldD ;
                D = DerivativePenalizer(s);
                D.penalize();
                dcsD1 = D.penalizedDerivative;
                
                

                v        = sum(sum(projectedFieldD))/(volfracD*elementNumberX*elementNumberY)-1;
                dv_dxs   = 1/(volfracD*elementNumberX*elementNumberY)*ones(elementNumberY,elementNumberX);

             %   Derivate the proyected field respected the filteredField
                s.beta = beta;
                s.filteredField = filteredField;

                s.eta =etaE;
                E = ProjectedFieldDerivator(s);
                E.compute();
                dxsE = E.derivatedProjectedField;

                s.eta =etaI;
                I = ProjectedFieldDerivator(s);
                I.compute();
                dxsI = I.derivatedProjectedField;

                s.eta =etaD;
                D = ProjectedFieldDerivator(s);
                D.compute();
                dxsD = D.derivatedProjectedField;

                
                %Derivate the filtered field by the field
                s.H = H;
                s.Hs = Hs;
                B = FilteredFieldDerivator(s);
                B.compute();
                derivedFilteredField = B.derivedFilteredField;
              
                % Calculate the cost derivated by the field
                s.derivedFilteredField = derivedFilteredField;

                s.derivedProyectedCost = dcsE1;
                s.derivatedProjectedField = dxsE;                
                E = CostDerivator(s);
                E.compute();
                dcsE1 = E.derivedCost;

                s.derivedProyectedCost = dcsI1;
                s.derivatedProjectedField = dxsI;                
                I = CostDerivator(s);
                I.compute();
                dcsI1 = I.derivedCost;

                s.derivedProyectedCost = dcsD1;
                s.derivatedProjectedField = dxsD;                
                D = CostDerivator(s);
                D.compute();
                dcsD1 = D.derivedCost;




                dvs_dxs(:)  = H*(dv_dxs(:).*dxsD(:)./Hs);

                % Compute the solver
                sC = [];
                c = CostAndConstraintArturo(sC);
                [f0val,df0dx,df0dx2] = c.computeCostValueAndGradient(elementNumberX,elementNumberY);
                [fval,dfdx,dfdx2] =  c.computeConstraintValueAndGradient(costE,costI,costD,v,dcsE1,dcsI1,dcsD1,dvs_dxs,cost);


                xval        = field(:);

                [xmma,~,zmma,~,~,~,~,~,~,low,upp] = mmasub2Arturo(numberRestriction,variableNumber,iter,xval,minDensity,maxDensity,xold1,xold2, ...
                    f0val,df0dx,df0dx2,fval,dfdx,dfdx2,low,upp,a0,mmaParameter.a,mmaParameter.e,mmaParameter.d);


                xold2  = xold1;
                xold1  = xval;


                field     = reshape(xmma,elementNumberY,elementNumberX);

                costChange = norm(xmma-xold,inf);

                filteredField(:) = (H*field(:))./Hs;

                %Project the new filtered field
             
                s.beta = beta;
                s.eta =etaE;
                s.filteredField =filteredField;
                E = FieldProjector(s);
                E.compute();
                projectedFieldE = E.projectedField;
    
                s.beta = beta;
                s.eta =etaI;
                s.filteredField =filteredField;
                I = FieldProjector(s);
                I.compute();
                projectedFieldI = I.projectedField;
    
                s.beta = beta;
                s.eta =etaD;
                s.filteredField =filteredField;
                D = FieldProjector(s);
                D.compute();
                projectedFieldD = D.projectedField;






                subplot(3,1,1);
                imagesc(-projectedFieldE); colormap(gray); caxis([-1 0]); axis off; axis equal; axis tight;
                subplot(3,1,2);
                imagesc(-projectedFieldI); colormap(gray); caxis([-1 0]); axis off; axis equal; axis tight;
                subplot(3,1,3);
                imagesc(-projectedFieldD); colormap(gray); caxis([-1 0]); axis off; axis equal; axis tight; pause(0.1)

                disp([' It.: ' sprintf('%4i',iter) ' Obj.: ' sprintf('%6.4f',zmma) ...
                    ' ui: '  sprintf('%12f', [costE costI costD])...
                    ' V: '   sprintf('%6.3f',sum(projectedFieldI(:))/(elementNumberX*elementNumberY)) ...
                    ' ch.: ' sprintf('%6.3f', costChange)])

                if (beta < 32) && ((iterbeta >= 50) || (costChange <= 0.01))
                    beta     = 2*beta;
                    iterbeta = 0;
                    costChange   = 1;
                end

                if (itervol >= 20 )
                    itervol = 0;
                    volfracD = volumenFraction*sum(projectedFieldD(:))/sum(projectedFieldI(:));
                end
                if iter == obj.iterations               
                    finish = true;
                end
                obj.projectedField.E = projectedFieldE;
                obj.projectedField.I = projectedFieldI;
                obj.projectedField.D = projectedFieldD;

            end

        end
    end
      

end