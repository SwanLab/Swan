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

%             %Matriz de rigidez (para elementos cuadrados):
%             A11     = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
%             A12     = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
%             B11     = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
%             B12     = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
%             elementalStiffnessMatrix      = t/(1-poissonCoefficient^2)/24*([A11 A12;A12' A11]+poissonCoefficient*[B11 B12;B12' B11]);
            
            s.elementType = 'square';
            s.t = t;
            s.poissonCoefficient = poissonCoefficient;
            B = ElementalStiffnessMatricesComputer(s);
            B.compute();
            elementalStiffnessMatrix = B.elementalStiffnessMatrix;


            s.elementNumberX =  elementNumberX;
            s.elementNumberY =  elementNumberY;
            B = GeometryComputer(s);
            B.compute();
            allDegress   = B.degress.all;
            freeDegress  = B.degress.free;
            fixedDegress = B.degress.fixed;
            conectivityMatrixMat = B.conectivityMatrixMat;


            
            iK      = reshape(kron(conectivityMatrixMat,ones(8,1))',64*elementNumberX*elementNumberY,1);
            jK      = reshape(kron(conectivityMatrixMat,ones(1,8))',64*elementNumberX*elementNumberY,1);
            output    = 2;

            %Inicializaciones: %Es una clase de desplazamientos; 2 dof por cada
            %elemento
            displacementE = zeros(2*(elementNumberY+1)*(elementNumberX+1),1);
            displacementI = zeros(2*(elementNumberY+1)*(elementNumberX+1),1);
            displacementD = zeros(2*(elementNumberY+1)*(elementNumberX+1),1);
            filteredField        = volumenFraction*ones(elementNumberY,elementNumberX);

            %Parámetros del MMA:
            % rowsolid  = 46:55;
            % columnsolid = 146:155;
            rowsolid = [];
            columnsolid = [];
            filteredField(rowsolid,columnsolid) = 1;


            minDensity    = zeros(elementNumberY,elementNumberX);
            maxDensity    =  ones(elementNumberY,elementNumberX);
            minDensity(rowsolid,columnsolid) = 0.99;


            minDensity = minDensity(:);
            maxDensity = maxDensity(:);

            % xmin    = zeros(nelx*nely,1);
            % xmax    =  ones(nelx*nely,1);
            costChange  = 1;
            numberRestriction       = 4;
            variableNumber       = elementNumberX*elementNumberY;
            low     = minDensity;
            upp     = maxDensity;
            a0      = 1;
            mmaParameter.a       = [1 1 1 0]';
            mmaParameter.e       = 1000*ones(numberRestriction,1);
            mmaParameter.d       = 0*ones(numberRestriction,1);
            xold1   = filteredField(:);
            xold2   = filteredField(:);

            %Preparación del filtro y de la proyección:
            iH = ones(elementNumberX*elementNumberY*(2*(ceil(minimunInfluenceRadios)-1)+1)^2,1);
            jH = ones(size(iH));
            sH = zeros(size(iH));
            k  = 0;
            for i1 = 1:elementNumberX
                for j1 = 1:elementNumberY
                    e1 = (i1-1)*elementNumberY+j1;
                    for i2 = max(i1-(ceil(minimunInfluenceRadios)-1),1):min(i1+(ceil(minimunInfluenceRadios)-1),elementNumberX)
                        for j2 = max(j1-(ceil(minimunInfluenceRadios)-1),1):min(j1+(ceil(minimunInfluenceRadios)-1),elementNumberY)
                            e2 = (i2-1)*elementNumberY+j2;
                            k = k+1;
                            iH(k) = e1;
                            jH(k) = e2;
                            sH(k) = max(0,minimunInfluenceRadios-sqrt((i1-i2)^2+(j1-j2)^2));
                        end
                    end
                end
            end
            H        = sparse(iH,jH,sH);
            Hs       = sum(H,2);

            eta      = 0.25;

            etaE     = 1-eta;
            etaI     = 0.5;
            etaD     = eta;
            beta     = 1;

            %Phys--> Definen la proyeccion, minimun lenght scale. Eta--> Controlar el
            %ancho de la escala.
            xTilde   = filteredField;
            %             xPhysE   = (tanh(beta*etaE) + tanh(beta*(xTilde-etaE)))/(tanh(beta*etaE) + tanh(beta*(1-etaE)));
            %             xPhysI   = (tanh(beta*etaI) + tanh(beta*(xTilde-etaI)))/(tanh(beta*etaI) + tanh(beta*(1-etaI)));
            %             xPhysD   = (tanh(beta*etaD) + tanh(beta*(xTilde-etaD)))/(tanh(beta*etaD) + tanh(beta*(1-etaD)));
            s.beta = beta;
            s.eta =etaE;
            s.filteredField =filteredField;
            E = FieldProjector(s);
            E.compute();
            xPhysE = E.projectedField;

            s.beta = beta;
            s.eta =etaI;
            s.filteredField =filteredField;
            I = FieldProjector(s);
            I.compute();
            xPhysI = I.projectedField;

            s.beta = beta;
            s.eta =etaD;
            s.filteredField =filteredField;
            D = FieldProjector(s);
            D.compute();
            xPhysD = D.projectedField;

            volfracD = volumenFraction*sum(xPhysD(:))/sum(xPhysI(:));

            iterbeta = 0;
            iter     = 0;
            itervol  = 0;
            finish = false;


           
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
            s.projectedField = xPhysE;
            
            B = FEMcomputer(s);
            B.compute();
            cte = B.cost;
            

            %Optimización
            while (costChange > 0.001) && (iter < 300) && (finish == false)

                %Actualización de variables al comienzo del bucle:
                iter     = iter     + 1;
                iterbeta = iterbeta + 1;
                itervol  = itervol  + 1;

                xold = filteredField(:);


             %Resolución del problema de elementos finitos:
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
                
                s.projectedField = xPhysE;                      
                E = FEMcomputer(s);
                E.compute();
                costE = E.cost;
                displacementE = E.displacement;

                s.projectedField = xPhysI;                      
                I = FEMcomputer(s);
                I.compute();
                costI = I.cost;
                displacementI = I.displacement;

                s.projectedField = xPhysD;                      
                D = FEMcomputer(s);
                D.compute();
                costD = D.cost;
                displacementD = D.displacement;
                           
              %Cálculo de las derivadas respecto a la variable estructural:

                s.elementNumberX = elementNumberX;
                s.elementNumberY = elementNumberY;
                s.elementalStiffnessMatrix = elementalStiffnessMatrix;
                s.conectivityMatrixMat = conectivityMatrixMat; 

                s.displacement =displacementE;
                E = CostDerivator(s);
                E.compute();
                cE1a = E.derivedCost;
                dcsE1 = -penalization*(elasticModuleNeutral-elasticModuleMinimun)*xPhysE.^(penalization-1).*cE1a;

                s.displacement =displacementI;
                I = CostDerivator(s);
                I.compute();
                cI1a = I.derivedCost;
                dcsI1 = -penalization*(elasticModuleNeutral-elasticModuleMinimun)*xPhysI.^(penalization-1).*cI1a;

                s.displacement =displacementD;
                D = CostDerivator(s);
                D.compute();
                cD1a = D.derivedCost;
                dcsD1 = -penalization*(elasticModuleNeutral-elasticModuleMinimun)*xPhysD.^(penalization-1).*cD1a;
                
                
           


                v        = sum(sum(xPhysD))/(volfracD*elementNumberX*elementNumberY)-1;
                dv_dxs   = 1/(volfracD*elementNumberX*elementNumberY)*ones(elementNumberY,elementNumberX);

                dxsE     = beta*(1 - (tanh(beta*(xTilde-etaE))).^2)/(tanh(beta*etaE) + tanh(beta*(1-etaE)));
                dcsE1(:) = H*(dcsE1(:).*dxsE(:)./Hs);

                dxsI     = beta*(1 - (tanh(beta*(xTilde-etaI))).^2)/(tanh(beta*etaI) + tanh(beta*(1-etaI)));
                dcsI1(:) = H*(dcsI1(:).*dxsI(:)./Hs);

                dxsD     = beta*(1 - (tanh(beta*(xTilde-etaD))).^2)/(tanh(beta*etaD) + tanh(beta*(1-etaD)));
                dcsD1(:) = H*(dcsD1(:).*dxsD(:)./Hs);

                dxsD     = beta*(1 - (tanh(beta*(xTilde-etaD))).^2)/(tanh(beta*etaD) + tanh(beta*(1-etaD)));

                dvs_dxs(:)  = H*(dv_dxs(:).*dxsD(:)./Hs);

                %MMA
                sC = [];
                c = CostAndConstraintArturo(sC);
                [f0val,df0dx,df0dx2] = c.computeCostValueAndGradient(elementNumberX,elementNumberY);
                [fval,dfdx,dfdx2] =  c.computeConstraintValueAndGradient(costE,costI,costD,v,dcsE1,dcsI1,dcsD1,dvs_dxs,cte);


                xval        = filteredField(:);

                %     [xmma,~,zmma,~,~,~,~,~,~,low,upp] = mmasub(m,n,iter,xval,xmin,xmax,xold1,xold2, ...
                %         f0val,df0dx,fval,dfdx,low,upp,a0,a,e,d);

                [xmma,~,zmma,~,~,~,~,~,~,low,upp] = mmasub2Arturo(numberRestriction,variableNumber,iter,xval,minDensity,maxDensity,xold1,xold2, ...
                    f0val,df0dx,df0dx2,fval,dfdx,dfdx2,low,upp,a0,mmaParameter.a,mmaParameter.e,mmaParameter.d);


                xold2  = xold1;
                xold1  = xval;


                filteredField     = reshape(xmma,elementNumberY,elementNumberX);

                costChange = norm(xmma-xold,inf);

                xTilde(:) = (H*filteredField(:))./Hs;
                xPhysE    = (tanh(beta*etaE) + tanh(beta*(xTilde-etaE)))/(tanh(beta*etaE) + tanh(beta*(1-etaE)));
                xPhysI    = (tanh(beta*etaI) + tanh(beta*(xTilde-etaI)))/(tanh(beta*etaI) + tanh(beta*(1-etaI)));
                xPhysD    = (tanh(beta*etaD) + tanh(beta*(xTilde-etaD)))/(tanh(beta*etaD) + tanh(beta*(1-etaD)));

                subplot(3,1,1);
                imagesc(-xPhysE); colormap(gray); caxis([-1 0]); axis off; axis equal; axis tight;
                subplot(3,1,2);
                imagesc(-xPhysI); colormap(gray); caxis([-1 0]); axis off; axis equal; axis tight;
                subplot(3,1,3);
                imagesc(-xPhysD); colormap(gray); caxis([-1 0]); axis off; axis equal; axis tight; pause(0.1)

                disp([' It.: ' sprintf('%4i',iter) ' Obj.: ' sprintf('%6.4f',zmma) ...
                    ' ui: '  sprintf('%12f', [costE costI costD])...
                    ' V: '   sprintf('%6.3f',sum(xPhysI(:))/(elementNumberX*elementNumberY)) ...
                    ' ch.: ' sprintf('%6.3f', costChange)])

                if (beta < 32) && ((iterbeta >= 50) || (costChange <= 0.01))
                    beta     = 2*beta;
                    iterbeta = 0;
                    costChange   = 1;
                end

                if (itervol >= 20 )
                    itervol = 0;
                    volfracD = volumenFraction*sum(xPhysD(:))/sum(xPhysI(:));
                end
                if iter == obj.iterations               
                    finish = true;
                end
                obj.projectedField.E = xPhysE;
                obj.projectedField.I = xPhysI;
                obj.projectedField.D = xPhysD;

            end

        end
    end
      

end