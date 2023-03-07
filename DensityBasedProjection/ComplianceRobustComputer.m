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


%Calculamos una normalización para el coste:
%             sKiI  = reshape(elementalStiffnessMatrix(:)*(elasticModuleMinimun + xPhysE(:)'.^penalization*(elasticModuleNeutral-elasticModuleMinimun)),64*elementNumberX*elementNumberY,1);
%             KiI   = sparse(iK,jK,sKiI);
%             KiI = (KiI+KiI')/2;
%             s.elementType = 'square';
%             s.t=t;
%             s.poissonCoefficient=poissonCoefficient;
%             s.elasticModuleMinimun=elasticModuleMinimun;
%             s.elasticModuleNeutral=elasticModuleNeutral;       
%             s.penalization=penalization;
%             s.elementNumberX = elementNumberX;
%             s.elementNumberY = elementNumberY;
%             s.projectedField= xPhysE;
%             s.conectivityMatrixMat = conectivityMatrixMat;
%             B = StifnessMatrixComputer(s);
%             B.compute();
%             KiI = B.globalStifnessMatrix;
% 
%             s.neumanCondition = neumanCondition;
%             s.output = output;
%             s.elementNumberX = elementNumberX;
%             s.elementNumberY = elementNumberY;
%             B = ForceComputer(s);
%             B.compute();
%             F = B.force;
% 

%             s.force = F;
%             s.globalStifnessMatrix =KiI;
%             s.elementNumberX =elementNumberX;
%             s.elementNumberY =elementNumberY;
%             s.freeDegress =freeDegress;
%             B = DisplacementComputer(s); 
%             B.compute()
%             UC = B.displacement;
%            

%             s.force = F;
%             s.displacement =UC;
%             B = CostComputer(s);
%             B.compute();
%             cte = B.cost;


           
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
                
                
           


                %     dcsE1(1:10,1:10)/cte
                %
                %     if iter == 2
                %     break;
                %     end

                %Derivadas numéricas
                %     dcnum = zeros(nely,nelx);
                %
                %     h = 1e-5;
                %
                %     costenum = zeros(nely,nelx);
                %
                %     for elx = 1:nelx
                %         for ely = 1:nely
                %             Um           = zeros(2*(nelx+1)*(nely+1),1);
                %             xsm          = xPhysE;
                %             xsm(ely,elx) = xsm(ely,elx) - h;
                %
                %             sKm  = reshape(ke(:)*(Emin + xsm(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
                %             Kim   = sparse(iK,jK,sKm); Kim = (Kim+Kim')/2;
                %
                %             Um(freedofs) = Kim(freedofs,freedofs)\F(freedofs);
                %
                %             costenum(ely,elx) = F'*Um;
                %         end
                %     end
                %
                %     dcnum = ((costE-costenum)/h)/cte;
                %     dcnum(1:10,1:10)
                %
                %     dcnum./dcsE1*cte

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
                f0val       = 0;
                df0dx       = zeros(elementNumberX*elementNumberY,1);
                df0dx2      = 0*df0dx;
                fval        = [costE; costI; costD; v];
                fval(1:3)   = fval(1:3)/cte;
                dfdx        = [ dcsE1(:)'; dcsI1(:)'; dcsD1(:)'; dvs_dxs(:)'];
                dfdx(1:3,:) = dfdx(1:3,:)/cte;
                dfdx2       = 0*dfdx;
                xval        = filteredField(:);

                %     [xmma,~,zmma,~,~,~,~,~,~,low,upp] = mmasub(m,n,iter,xval,xmin,xmax,xold1,xold2, ...
                %         f0val,df0dx,fval,dfdx,low,upp,a0,a,e,d);

                [xmma,~,zmma,~,~,~,~,~,~,low,upp] = ComplianceRobustComputer.mmasub2(numberRestriction,variableNumber,iter,xval,minDensity,maxDensity,xold1,xold2, ...
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
    methods (Static, Access = private)
        function [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = ...
                mmasub2(m,n,iter,xval,xmin,xmax,xold1,xold2, ...
                f0val,df0dx,df0dx2,fval,dfdx,dfdx2,low,upp,a0,a,c,d);
            %
            %    Written in May 1999 by
            %    Krister Svanberg <krille@math.kth.se>
            %    Department of Mathematics
            %    SE-10044 Stockholm, Sweden.
            %
            %    This function mmasub performs one MMA-iteration, aimed at
            %    solving the nonlinear programming problem:
            %
            %      Minimize  f_0(x) + a_0*z + sum( c_i*y_i + 0.5*d_i*(y_i)^2 )
            %    subject to  f_i(x) - a_i*z - y_i <= 0,  i = 1,...,m
            %                xmin_j <= x_j <= xmax_j,    j = 1,...,n
            %                z >= 0,   y_i >= 0,         i = 1,...,m
            %*** INPUT:
            %
            %   m    = The number of general constraints.
            %   n    = The number of variables x_j.
            %  iter  = Current iteration number ( =1 the first time mmasub is called).
            %  xval  = Column vector with the current values of the variables x_j.
            %  xmin  = Column vector with the lower bounds for the variables x_j.
            %  xmax  = Column vector with the upper bounds for the variables x_j.
            %  xold1 = xval, one iteration ago (provided that iter>1).
            %  xold2 = xval, two iterations ago (provided that iter>2).
            %  f0val = The value of the objective function f_0 at xval.
            %  df0dx = Column vector with the derivatives of the objective function
            %          f_0 with respect to the variables x_j, calculated at xval.
            % df0dx2 = Column vector with the non-mixed second derivatives of the
            %          objective function f_0 with respect to the variables x_j,
            %          calculated at xval. df0dx2(j) = the second derivative
            %          of f_0 with respect to x_j (twice).
            %          Important note: If second derivatives are not available,
            %          simply let df0dx2 = 0*df0dx.
            %  fval  = Column vector with the values of the constraint functions f_i,
            %          calculated at xval.
            %  dfdx  = (m x n)-matrix with the derivatives of the constraint functions
            %          f_i with respect to the variables x_j, calculated at xval.
            %          dfdx(i,j) = the derivative of f_i with respect to x_j.
            %  dfdx2 = (m x n)-matrix with the non-mixed second derivatives of the
            %          constraint functions f_i with respect to the variables x_j,
            %          calculated at xval. dfdx2(i,j) = the second derivative
            %          of f_i with respect to x_j (twice).
            %          Important note: If second derivatives are not available,
            %          simply let dfdx2 = 0*dfdx.
            %  low   = Column vector with the lower asymptotes from the previous
            %          iteration (provided that iter>1).
            %  upp   = Column vector with the upper asymptotes from the previous
            %          iteration (provided that iter>1).
            %  a0    = The constants a_0 in the term a_0*z.
            %  a     = Column vector with the constants a_i in the terms a_i*z.
            %  c     = Column vector with the constants c_i in the terms c_i*y_i.
            %  d     = Column vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.
            %
            %*** OUTPUT:
            %
            %  xmma  = Column vector with the optimal values of the variables x_j
            %          in the current MMA subproblem.
            %  ymma  = Column vector with the optimal values of the variables y_i
            %          in the current MMA subproblem.
            %  zmma  = Scalar with the optimal value of the variable z
            %          in the current MMA subproblem.
            %  lam   = Lagrange multipliers for the m general MMA constraints.
            %  xsi   = Lagrange multipliers for the n constraints alfa_j - x_j <= 0.
            %  eta   = Lagrange multipliers for the n constraints x_j - beta_j <= 0.
            %   mu   = Lagrange multipliers for the m constraints -y_i <= 0.
            %  zet   = Lagrange multiplier for the single constraint -z <= 0.
            %   s    = Slack variables for the m general MMA constraints.
            %  low   = Column vector with the lower asymptotes, calculated and used
            %          in the current MMA subproblem.
            %  upp   = Column vector with the upper asymptotes, calculated and used
            %          in the current MMA subproblem.
            %
            epsimin = sqrt(m+n)*10^(-9);
            feps = 0.000001;
            asyinit = 0.2;
            asyincr = 1.1;
            asydecr = 0.65;
            albefa = 0.1;
            een = ones(n,1);
            zeron = zeros(n,1);

            % Calculation of the asymptotes low and upp :
            if iter < 2.5
                low = xval - asyinit*(xmax-xmin);
                upp = xval + asyinit*(xmax-xmin);
            else
                zzz = (xval-xold1).*(xold1-xold2);
                factor = een;
                factor(zzz > 0) = asyincr;
                factor(zzz < 0) = asydecr;
                low = xval - factor.*(xold1 - low);
                upp = xval + factor.*(upp - xold1);
            end

            % Calculation of the bounds alfa and beta :
            zzz = low + albefa*(xval-low);
            alfa = max(zzz,xmin);
            zzz = upp - albefa*(upp-xval);
            beta = min(zzz,xmax);

            % Calculations of p0, q0, P, Q and b.


            ux1 = upp-xval;
            ux2 = ux1.*ux1;
            ux3 = ux2.*ux1;
            xl1 = xval-low;
            xl2 = xl1.*xl1;
            xl3 = xl2.*xl1;
            ul1 = upp-low;
            ulinv1 = een./ul1;
            uxinv1 = een./ux1;
            xlinv1 = een./xl1;
            uxinv3 = een./ux3;
            xlinv3 = een./xl3;
            diap = (ux3.*xl1)./(2*ul1);
            diaq = (ux1.*xl3)./(2*ul1);
            p0 = zeron;
            p0(df0dx > 0) = df0dx(df0dx > 0);
            p0 = p0 + 0.001*abs(df0dx) + feps*ulinv1;
            p0 = p0.*ux2;
            q0 = zeron;
            q0(df0dx < 0) = -df0dx(df0dx < 0);
            q0 = q0 + 0.001*abs(df0dx) + feps*ulinv1;
            q0 = q0.*xl2;
            dg0dx2 = 2*(p0./ux3 + q0./xl3);
            del0 = df0dx2 - dg0dx2;
            delpos0 = zeron;
            delpos0(del0 > 0) = del0(del0 > 0);
            p0 = p0 + delpos0.*diap;
            q0 = q0 + delpos0.*diaq;
            P = zeros(m,n);
            P(dfdx > 0) = dfdx(dfdx > 0);
            P = P * diag(ux2);
            Q = zeros(m,n);
            Q(dfdx < 0) = -dfdx(dfdx < 0);
            Q = Q * diag(xl2);
            dgdx2 = 2*(P*diag(uxinv3) + Q*diag(xlinv3));
            del = dfdx2 - dgdx2;
            delpos = zeros(m,n);
            delpos(del > 0) = del(del > 0);
            P = P + delpos*diag(diap);
            Q = Q + delpos*diag(diaq);
            b = P*uxinv1 + Q*xlinv1 - fval ;


            %%% Solving the subproblem by a primal-dual Newton method
            [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s] = ...
                ComplianceRobustComputer.subsolv2(m,n,epsimin,low,upp,alfa,beta,p0,q0,P,Q,a0,a,b,c,d);

        end


        function [xmma,ymma,zmma,lamma,xsimma,etamma,mumma,zetmma,smma] = ...
                subsolv2(m,n,epsimin,low,upp,alfa,beta,p0,q0,P,Q,a0,a,b,c,d)
            %
            %    Written in May 1999 by
            %    Krister Svanberg <krille@math.kth.se>
            %    Department of Mathematics
            %    SE-10044 Stockholm, Sweden.
            %
            % This function subsolv solves the MMA subproblem:
            %
            % minimize   SUM[ p0j/(uppj-xj) + q0j/(xj-lowj) ] + a0*z +
            %          + SUM[ ci*yi + 0.5*di*(yi)^2 ],
            %
            % subject to SUM[ pij/(uppj-xj) + qij/(xj-lowj) ] - ai*z - yi <= bi,
            %            alfaj <=  xj <=  betaj,  yi >= 0,  z >= 0.
            %
            % Input:  m, n, low, upp, alfa, beta, p0, q0, P, Q, a0, a, b, c, d.
            % Output: xmma,ymma,zmma, slack variables and Lagrange multiplers.
            %
            een = ones(n,1);
            eem = ones(m,1);
            epsi = 1;
            epsvecn = epsi*een;
            epsvecm = epsi*eem;
            x = 0.5*(alfa+beta);
            y = eem;
            z = 1;
            lam = eem;
            xsi = een./(x-alfa);
            xsi = max(xsi,een);
            eta = een./(beta-x);
            eta = max(eta,een);
            mu  = max(eem,0.5*c);
            zet = 1;
            s = eem;
            itera = 0;

            while epsi > epsimin
                epsvecn = epsi*een;
                epsvecm = epsi*eem;
                ux1 = upp-x;
                xl1 = x-low;
                ux2 = ux1.*ux1;
                xl2 = xl1.*xl1;
                uxinv1 = een./ux1;
                xlinv1 = een./xl1;

                plam = p0 + P'*lam ;
                qlam = q0 + Q'*lam ;
                gvec = P*uxinv1 + Q*xlinv1;
                dpsidx = plam./ux2 - qlam./xl2 ;

                rex = dpsidx - xsi + eta;
                rey = c + d.*y - mu - lam;
                rez = a0 - zet - a'*lam;
                relam = gvec - a*z - y + s - b;
                rexsi = xsi.*(x-alfa) - epsvecn;
                reeta = eta.*(beta-x) - epsvecn;
                remu = mu.*y - epsvecm;
                rezet = zet*z - epsi;
                res = lam.*s - epsvecm;

                residu1 = [rex' rey' rez]';
                residu2 = [relam' rexsi' reeta' remu' rezet res']';
                residu = [residu1' residu2']';
                residunorm = sqrt(residu'*residu);
                residumax = max(abs(residu));

                ittt = 0;
                while residumax > 0.9*epsi && ittt < 100
                    ittt=ittt + 1;
                    itera=itera + 1;

                    ux1 = upp-x;
                    xl1 = x-low;
                    ux2 = ux1.*ux1;
                    xl2 = xl1.*xl1;
                    ux3 = ux1.*ux2;
                    xl3 = xl1.*xl2;
                    uxinv1 = een./ux1;
                    xlinv1 = een./xl1;
                    uxinv2 = een./ux2;
                    xlinv2 = een./xl2;
                    plam = p0 + P'*lam ;
                    qlam = q0 + Q'*lam ;
                    gvec = P*uxinv1 + Q*xlinv1;
                    GG = P*spdiags(uxinv2,0,n,n) - Q*spdiags(xlinv2,0,n,n);
                    dpsidx = plam./ux2 - qlam./xl2 ;
                    delx = dpsidx - epsvecn./(x-alfa) + epsvecn./(beta-x);
                    dely = c + d.*y - lam - epsvecm./y;
                    delz = a0 - a'*lam - epsi/z;
                    dellam = gvec - a*z - y - b + epsvecm./lam;
                    diagx = plam./ux3 + qlam./xl3;
                    diagx = 2*diagx + xsi./(x-alfa) + eta./(beta-x);
                    diagxinv = een./diagx;
                    diagy = d + mu./y;
                    diagyinv = eem./diagy;
                    diaglam = s./lam;
                    diaglamyi = diaglam+diagyinv;

                    if m < n
                        blam = dellam + dely./diagy - GG*(delx./diagx);
                        bb = [blam' delz]';
                        Alam = spdiags(diaglamyi,0,m,m) + GG*spdiags(diagxinv,0,n,n)*GG';
                        AA = [Alam     a
                            a'    -zet/z ];
                        solut = AA\bb;
                        dlam = solut(1:m);
                        dz = solut(m+1);
                        dx = -delx./diagx - (GG'*dlam)./diagx;
                    else
                        diaglamyiinv = eem./diaglamyi;
                        dellamyi = dellam + dely./diagy;
                        Axx = spdiags(diagx,0,n,n) + GG'*spdiags(diaglamyiinv,0,m,m)*GG;
                        azz = zet/z + a'*(a./diaglamyi);
                        axz = -GG'*(a./diaglamyi);
                        bx = delx + GG'*(dellamyi./diaglamyi);
                        bz  = delz - a'*(dellamyi./diaglamyi);
                        AA = [Axx   axz
                            axz'  azz ];
                        bb = [-bx' -bz]';
                        solut = AA\bb;
                        dx  = solut(1:n);
                        dz = solut(n+1);
                        dlam = (GG*dx)./diaglamyi - dz*(a./diaglamyi) + dellamyi./diaglamyi;
                    end

                    dy = -dely./diagy + dlam./diagy;
                    dxsi = -xsi + epsvecn./(x-alfa) - (xsi.*dx)./(x-alfa);
                    deta = -eta + epsvecn./(beta-x) + (eta.*dx)./(beta-x);
                    dmu  = -mu + epsvecm./y - (mu.*dy)./y;
                    dzet = -zet + epsi/z - zet*dz/z;
                    ds   = -s + epsvecm./lam - (s.*dlam)./lam;
                    xx  = [ y'  z  lam'  xsi'  eta'  mu'  zet  s']';
                    dxx = [dy' dz dlam' dxsi' deta' dmu' dzet ds']';

                    stepxx = -1.01*dxx./xx;
                    stmxx  = max(stepxx);
                    stepalfa = -1.01*dx./(x-alfa);
                    stmalfa = max(stepalfa);
                    stepbeta = 1.01*dx./(beta-x);
                    stmbeta = max(stepbeta);
                    stmalbe  = max(stmalfa,stmbeta);
                    stmalbexx = max(stmalbe,stmxx);
                    stminv = max(stmalbexx,1);
                    steg = 1/stminv;

                    xold   =   x;
                    yold   =   y;
                    zold   =   z;
                    lamold =  lam;
                    xsiold =  xsi;
                    etaold =  eta;
                    muold  =  mu;
                    zetold =  zet;
                    sold   =   s;

                    itto = 0;
                    resinew = 2*residunorm;
                    while resinew > residunorm && itto < 50
                        itto = itto+1;

                        x   =   xold + steg*dx;
                        y   =   yold + steg*dy;
                        z   =   zold + steg*dz;
                        lam = lamold + steg*dlam;
                        xsi = xsiold + steg*dxsi;
                        eta = etaold + steg*deta;
                        mu  = muold  + steg*dmu;
                        zet = zetold + steg*dzet;
                        s   =   sold + steg*ds;
                        ux1 = upp-x;
                        xl1 = x-low;
                        ux2 = ux1.*ux1;
                        xl2 = xl1.*xl1;
                        uxinv1 = een./ux1;
                        xlinv1 = een./xl1;
                        plam = p0 + P'*lam ;
                        qlam = q0 + Q'*lam ;
                        gvec = P*uxinv1 + Q*xlinv1;
                        dpsidx = plam./ux2 - qlam./xl2 ;

                        rex = dpsidx - xsi + eta;
                        rey = c + d.*y - mu - lam;
                        rez = a0 - zet - a'*lam;
                        relam = gvec - a*z - y + s - b;
                        rexsi = xsi.*(x-alfa) - epsvecn;
                        reeta = eta.*(beta-x) - epsvecn;
                        remu = mu.*y - epsvecm;
                        rezet = zet*z - epsi;
                        res = lam.*s - epsvecm;

                        residu1 = [rex' rey' rez]';
                        residu2 = [relam' rexsi' reeta' remu' rezet res']';
                        residu = [residu1' residu2']';
                        resinew = sqrt(residu'*residu);
                        steg = steg/2;
                    end
                    residunorm=resinew;
                    residumax = max(abs(residu));
                    steg = 2*steg;
                end
                epsi = 0.1*epsi;
            end

            xmma   =   x;
            ymma   =   y;
            zmma   =   z;
            lamma =  lam;
            xsimma =  xsi;
            etamma =  eta;
            mumma  =  mu;
            zetmma =  zet;
            smma   =   s;

        end







        function [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = ...
                mmasub(m,n,iter,xval,xmin,xmax,xold1,xold2, ...
                ~,df0dx,fval,dfdx,low,upp,a0,a,c,d)
            %
            %    Version September 2007 (and a small change August 2008)
            %
            %    Krister Svanberg <krille@math.kth.se>
            %    Department of Mathematics, SE-10044 Stockholm, Sweden.
            %
            %    This function mmasub performs one MMA-iteration, aimed at
            %    solving the nonlinear programming problem:
            %
            %      Minimize  f_0(x) + a_0*z + sum( c_i*y_i + 0.5*d_i*(y_i)^2 )
            %    subject to  f_i(x) - a_i*z - y_i <= 0,  i = 1,...,m
            %                xmin_j <= x_j <= xmax_j,    j = 1,...,n
            %                z >= 0,   y_i >= 0,         i = 1,...,m
            %*** INPUT:
            %
            %   m    = The number of general constraints.
            %   n    = The number of variables x_j.
            %  iter  = Current iteration number ( =1 the first time mmasub is called).
            %  xval  = Column vector with the current values of the variables x_j.
            %  xmin  = Column vector with the lower bounds for the variables x_j.
            %  xmax  = Column vector with the upper bounds for the variables x_j.
            %  xold1 = xval, one iteration ago (provided that iter>1).
            %  xold2 = xval, two iterations ago (provided that iter>2).
            %  f0val = The value of the objective function f_0 at xval.
            %  df0dx = Column vector with the derivatives of the objective function
            %          f_0 with respect to the variables x_j, calculated at xval.
            %  fval  = Column vector with the values of the constraint functions f_i,
            %          calculated at xval.
            %  dfdx  = (m x n)-matrix with the derivatives of the constraint functions
            %          f_i with respect to the variables x_j, calculated at xval.
            %          dfdx(i,j) = the derivative of f_i with respect to x_j.
            %  low   = Column vector with the lower asymptotes from the previous
            %          iteration (provided that iter>1).
            %  upp   = Column vector with the upper asymptotes from the previous
            %          iteration (provided that iter>1).
            %  a0    = The constants a_0 in the term a_0*z.
            %  a     = Column vector with the constants a_i in the terms a_i*z.
            %  c     = Column vector with the constants c_i in the terms c_i*y_i.
            %  d     = Column vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.
            %
            %*** OUTPUT:
            %
            %  xmma  = Column vector with the optimal values of the variables x_j
            %          in the current MMA subproblem.
            %  ymma  = Column vector with the optimal values of the variables y_i
            %          in the current MMA subproblem.
            %  zmma  = Scalar with the optimal value of the variable z
            %          in the current MMA subproblem.
            %  lam   = Lagrange multipliers for the m general MMA constraints.
            %  xsi   = Lagrange multipliers for the n constraints alfa_j - x_j <= 0.
            %  eta   = Lagrange multipliers for the n constraints x_j - beta_j <= 0.
            %   mu   = Lagrange multipliers for the m constraints -y_i <= 0.
            %  zet   = Lagrange multiplier for the single constraint -z <= 0.
            %   s    = Slack variables for the m general MMA constraints.
            %  low   = Column vector with the lower asymptotes, calculated and used
            %          in the current MMA subproblem.
            %  upp   = Column vector with the upper asymptotes, calculated and used
            %          in the current MMA subproblem.
            %
            %epsimin = sqrt(m+n)*10^(-9);
            epsimin = 10^(-7);
            raa0 = 0.00001;
            move = .20;
            albefa = 0.1;
            asyinit = 0.5;
            asyincr = 1.2;
            asydecr = 0.7;
            eeen = ones(n,1);
            eeem = ones(m,1);
            zeron = zeros(n,1);

            % Calculation of the asymptotes low and upp :
            if iter < 2.5
                low = xval - asyinit*(xmax-xmin);
                upp = xval + asyinit*(xmax-xmin);
            else
                zzz = (xval-xold1).*(xold1-xold2);
                factor = eeen;
                factor(zzz > 0) = asyincr;
                factor(zzz < 0) = asydecr;
                low = xval - factor.*(xold1 - low);
                upp = xval + factor.*(upp - xold1);
                lowmin = xval - 10*(xmax-xmin);
                lowmax = xval - 0.01*(xmax-xmin);
                uppmin = xval + 0.01*(xmax-xmin);
                uppmax = xval + 10*(xmax-xmin);
                low = max(low,lowmin);
                low = min(low,lowmax);
                upp = min(upp,uppmax);
                upp = max(upp,uppmin);
            end

            % Calculation of the bounds alfa and beta :

            zzz1 = low + albefa*(xval-low);
            zzz2 = xval - move*(xmax-xmin);
            zzz  = max(zzz1,zzz2);
            alfa = max(zzz,xmin);
            zzz1 = upp - albefa*(upp-xval);
            zzz2 = xval + move*(xmax-xmin);
            zzz  = min(zzz1,zzz2);
            beta = min(zzz,xmax);

            % Calculations of p0, q0, P, Q and b.

            xmami = xmax-xmin;
            xmamieps = 0.00001*eeen;
            xmami = max(xmami,xmamieps);
            xmamiinv = eeen./xmami;
            ux1 = upp-xval;
            ux2 = ux1.*ux1;
            xl1 = xval-low;
            xl2 = xl1.*xl1;
            uxinv = eeen./ux1;
            xlinv = eeen./xl1;
            %
            p0 = zeron;
            q0 = zeron;
            p0 = max(df0dx,0);
            q0 = max(-df0dx,0);
            %p0(find(df0dx > 0)) = df0dx(find(df0dx > 0));
            %q0(find(df0dx < 0)) = -df0dx(find(df0dx < 0));
            pq0 = 0.001*(p0 + q0) + raa0*xmamiinv;
            p0 = p0 + pq0;
            q0 = q0 + pq0;
            p0 = p0.*ux2;
            q0 = q0.*xl2;
            %
            P = sparse(m,n);
            Q = sparse(m,n);
            P = max(dfdx,0);
            Q = max(-dfdx,0);
            %P(find(dfdx > 0)) = dfdx(find(dfdx > 0));
            %Q(find(dfdx < 0)) = -dfdx(find(dfdx < 0));
            PQ = 0.001*(P + Q) + raa0*eeem*xmamiinv';
            P = P + PQ;
            Q = Q + PQ;
            P = P * spdiags(ux2,0,n,n);
            Q = Q * spdiags(xl2,0,n,n);
            b = P*uxinv + Q*xlinv - fval ;
            %
            %%% Solving the subproblem by a primal-dual Newton method
            [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s] = ...
                ComplianceRobustComputer.subsolv(m,n,epsimin,low,upp,alfa,beta,p0,q0,P,Q,a0,a,b,c,d);

        end

        %-------------------------------------------------------------
        %    This is the file subsolv.m
        %
        %    Version Dec 2006.
        %    Krister Svanberg <krille@math.kth.se>
        %    Department of Mathematics, KTH,
        %    SE-10044 Stockholm, Sweden.
        %
        function [xmma,ymma,zmma,lamma,xsimma,etamma,mumma,zetmma,smma] = ...
                subsolv(m,n,epsimin,low,upp,alfa,beta,p0,q0,P,Q,a0,a,b,c,d)
            %
            % This function subsolv solves the MMA subproblem:
            %
            % minimize   SUM[ p0j/(uppj-xj) + q0j/(xj-lowj) ] + a0*z +
            %          + SUM[ ci*yi + 0.5*di*(yi)^2 ],
            %
            % subject to SUM[ pij/(uppj-xj) + qij/(xj-lowj) ] - ai*z - yi <= bi,
            %            alfaj <=  xj <=  betaj,  yi >= 0,  z >= 0.
            %
            % Input:  m, n, low, upp, alfa, beta, p0, q0, P, Q, a0, a, b, c, d.
            % Output: xmma,ymma,zmma, slack variables and Lagrange multiplers.
            %
            een = ones(n,1);
            eem = ones(m,1);
            epsi = 1;
            %epsvecn = epsi*een;
            %epsvecm = epsi*eem;
            x = 0.5*(alfa+beta);
            y = eem;
            z = 1;
            lam = eem;
            xsi = een./(x-alfa);
            xsi = max(xsi,een);
            eta = een./(beta-x);
            eta = max(eta,een);
            mu  = max(eem,0.5*c);
            zet = 1;
            s = eem;
            itera = 0;
            while epsi > epsimin
                epsvecn = epsi*een;
                epsvecm = epsi*eem;
                ux1 = upp-x;
                xl1 = x-low;
                ux2 = ux1.*ux1;
                xl2 = xl1.*xl1;
                uxinv1 = een./ux1;
                xlinv1 = een./xl1;
                plam = p0 + P'*lam ;
                qlam = q0 + Q'*lam ;
                gvec = P*uxinv1 + Q*xlinv1;
                dpsidx = plam./ux2 - qlam./xl2 ;
                rex = dpsidx - xsi + eta;
                rey = c + d.*y - mu - lam;
                rez = a0 - zet - a'*lam;
                relam = gvec - a*z - y + s - b;
                rexsi = xsi.*(x-alfa) - epsvecn;
                reeta = eta.*(beta-x) - epsvecn;
                remu = mu.*y - epsvecm;
                rezet = zet*z - epsi;
                res = lam.*s - epsvecm;
                residu1 = [rex' rey' rez]';
                residu2 = [relam' rexsi' reeta' remu' rezet res']';
                residu = [residu1' residu2']';
                residunorm = sqrt(residu'*residu);
                residumax = max(abs(residu));
                ittt = 0;
                while residumax > 0.9*epsi && ittt < 200
                    ittt=ittt + 1;
                    itera=itera + 1;
                    ux1 = upp-x;
                    xl1 = x-low;
                    ux2 = ux1.*ux1;
                    xl2 = xl1.*xl1;
                    ux3 = ux1.*ux2;
                    xl3 = xl1.*xl2;
                    uxinv1 = een./ux1;
                    xlinv1 = een./xl1;
                    uxinv2 = een./ux2;
                    xlinv2 = een./xl2;
                    plam = p0 + P'*lam ;
                    qlam = q0 + Q'*lam ;
                    gvec = P*uxinv1 + Q*xlinv1;
                    GG = P*spdiags(uxinv2,0,n,n) - Q*spdiags(xlinv2,0,n,n);
                    dpsidx = plam./ux2 - qlam./xl2 ;
                    delx = dpsidx - epsvecn./(x-alfa) + epsvecn./(beta-x);
                    dely = c + d.*y - lam - epsvecm./y;
                    delz = a0 - a'*lam - epsi/z;
                    dellam = gvec - a*z - y - b + epsvecm./lam;
                    diagx = plam./ux3 + qlam./xl3;
                    diagx = 2*diagx + xsi./(x-alfa) + eta./(beta-x);
                    diagxinv = een./diagx;
                    diagy = d + mu./y;
                    diagyinv = eem./diagy;
                    diaglam = s./lam;
                    diaglamyi = diaglam+diagyinv;
                    if m < n
                        blam = dellam + dely./diagy - GG*(delx./diagx);
                        bb = [blam' delz]';
                        Alam = spdiags(diaglamyi,0,m,m) + GG*spdiags(diagxinv,0,n,n)*GG';
                        AA = [Alam     a
                            a'    -zet/z ];
                        solut = AA\bb;
                        dlam = solut(1:m);
                        dz = solut(m+1);
                        dx = -delx./diagx - (GG'*dlam)./diagx;
                    else
                        diaglamyiinv = eem./diaglamyi;
                        dellamyi = dellam + dely./diagy;
                        Axx = spdiags(diagx,0,n,n) + GG'*spdiags(diaglamyiinv,0,m,m)*GG;
                        azz = zet/z + a'*(a./diaglamyi);
                        axz = -GG'*(a./diaglamyi);
                        bx = delx + GG'*(dellamyi./diaglamyi);
                        bz  = delz - a'*(dellamyi./diaglamyi);
                        AA = [Axx   axz
                            axz'  azz ];
                        bb = [-bx' -bz]';
                        solut = AA\bb;
                        dx  = solut(1:n);
                        dz = solut(n+1);
                        dlam = (GG*dx)./diaglamyi - dz*(a./diaglamyi) + dellamyi./diaglamyi;
                    end
                    %
                    dy = -dely./diagy + dlam./diagy;
                    dxsi = -xsi + epsvecn./(x-alfa) - (xsi.*dx)./(x-alfa);
                    deta = -eta + epsvecn./(beta-x) + (eta.*dx)./(beta-x);
                    dmu  = -mu + epsvecm./y - (mu.*dy)./y;
                    dzet = -zet + epsi/z - zet*dz/z;
                    ds   = -s + epsvecm./lam - (s.*dlam)./lam;
                    xx  = [ y'  z  lam'  xsi'  eta'  mu'  zet  s']';
                    dxx = [dy' dz dlam' dxsi' deta' dmu' dzet ds']';
                    %
                    stepxx = -1.01*dxx./xx;
                    stmxx  = max(stepxx);
                    stepalfa = -1.01*dx./(x-alfa);
                    stmalfa = max(stepalfa);
                    stepbeta = 1.01*dx./(beta-x);
                    stmbeta = max(stepbeta);
                    stmalbe  = max(stmalfa,stmbeta);
                    stmalbexx = max(stmalbe,stmxx);
                    stminv = max(stmalbexx,1);
                    steg = 1/stminv;
                    %
                    xold   =   x;
                    yold   =   y;
                    zold   =   z;
                    lamold =  lam;
                    xsiold =  xsi;
                    etaold =  eta;
                    muold  =  mu;
                    zetold =  zet;
                    sold   =   s;
                    %
                    itto = 0;
                    resinew = 2*residunorm;
                    while resinew > residunorm && itto < 50
                        itto = itto+1;
                        x   =   xold + steg*dx;
                        y   =   yold + steg*dy;
                        z   =   zold + steg*dz;
                        lam = lamold + steg*dlam;
                        xsi = xsiold + steg*dxsi;
                        eta = etaold + steg*deta;
                        mu  = muold  + steg*dmu;
                        zet = zetold + steg*dzet;
                        s   =   sold + steg*ds;
                        ux1 = upp-x;
                        xl1 = x-low;
                        ux2 = ux1.*ux1;
                        xl2 = xl1.*xl1;
                        uxinv1 = een./ux1;
                        xlinv1 = een./xl1;
                        plam = p0 + P'*lam ;
                        qlam = q0 + Q'*lam ;
                        gvec = P*uxinv1 + Q*xlinv1;
                        dpsidx = plam./ux2 - qlam./xl2 ;
                        rex = dpsidx - xsi + eta;
                        rey = c + d.*y - mu - lam;
                        rez = a0 - zet - a'*lam;
                        relam = gvec - a*z - y + s - b;
                        rexsi = xsi.*(x-alfa) - epsvecn;
                        reeta = eta.*(beta-x) - epsvecn;
                        remu = mu.*y - epsvecm;
                        rezet = zet*z - epsi;
                        res = lam.*s - epsvecm;
                        residu1 = [rex' rey' rez]';
                        residu2 = [relam' rexsi' reeta' remu' rezet res']';
                        residu = [residu1' residu2']';
                        resinew = sqrt(residu'*residu);
                        steg = steg/2;
                    end
                    residunorm=resinew;
                    residumax = max(abs(residu));
                    steg = 2*steg;
                end
                %if ittt > 198
                %  epsi;
                %  ittt;
                %end
                epsi = 0.1*epsi;
            end
            xmma   =   x;
            ymma   =   y;
            zmma   =   z;
            lamma =  lam;
            xsimma =  xsi;
            etamma =  eta;
            mumma  =  mu;
            zetmma =  zet;
            smma   =   s;
            %-------------------------------------------------------------
        end
    end

end