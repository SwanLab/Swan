classdef TopologicalDerivativeComputer < handle

    properties (Access = public)
        dt
        dC
        strain
        TD
        tgamma
    end

    properties (Access = private)
        p
        tensor
        t
        mesh
        psi
        designVariable
        penalization
        penalty
        volfrac
        auglag
        max_vol
        energy0
        voltarget
        nMat
        mat
        gamma
        tE
        beta
        alpha
        la0 
        mu0
        young
        stress
        U
        volume
    end

    methods (Access = public)

        function obj = TopologicalDerivativeComputer(cParams)
            obj.init(cParams)
            obj.computeGamma();
            obj.computeTgamma();
            obj.computeBetaAndAlpha();
            obj.computeLameParameters();
            obj.computeStressAndStrain();
            obj.computeTopologicalDerivative();
            %obj.computeVolumeConstraintinDT();
            %obj.smoothTopologicalDerivative();
        end
    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.p = cParams.meshSeba.p; 
            obj.t = cParams.meshSeba.t;
            obj.mesh = cParams.mesh;
            obj.penalization = cParams.penalization;
            obj.penalty = cParams.penalty;
            obj.volfrac = cParams.volfrac;
            obj.auglag = cParams.auglag;
            obj.max_vol = cParams.max_vol; 
            obj.energy0 = cParams.energy0;
            obj.voltarget = obj.max_vol*obj.volfrac;
            obj.nMat = cParams.nMat;
            obj.mat{1} = cParams.mat.A;
            obj.mat{2} = cParams.mat.B;
            obj.mat{3} = cParams.mat.C;
            obj.mat{4} = cParams.mat.D;
            obj.psi = cParams.psi;   
            obj.designVariable = cParams.designVariable;
            obj.U = cParams.U;
            obj.volume = cParams.volume;
        
        end

        function computeGamma(obj)
            E1 = obj.mat{1}.young;
            E2 = obj.mat{2}.young;
            E3 = obj.mat{3}.young;
            E4 = obj.mat{4}.young;

            obj.young = [E1 E2 E3 E4];

            obj.gamma(1) = E1/E1;
            obj.gamma(2) = E2/E1;
            obj.gamma(3) = E3/E1;
            obj.gamma(4) = E4/E1;    
        end

        function computeTgamma(obj)
            s.psi = obj.psi;
            s.p = obj.p;
            s.t = obj.t; 
            s.designVariable = obj.designVariable;
            s.m = obj.mesh;
            
            charfun = CharacteristicFunctionComputer(s); 
            [~,tfi] = charfun.computeFiandTfi();

            obj.tgamma = obj.gamma*tfi; %Mixed formulation method
            E1 = obj.mat{1}.young;
            obj.tE = E1*obj.tgamma; 
        end

        function computeBetaAndAlpha(obj)
             
            nu = obj.mat{1}.nu;
            obj.beta = (1+nu)/(1-nu); 
            obj.alpha = (3-nu)/(1+nu);
        end

        function computeLameParameters(obj)
            nu = obj.mat{1}.nu;
            E1 = obj.mat{1}.young;
        
            obj.la0 = nu.*E1./((1+nu).*(1-2.*nu)); obj.mu0 = E1./(2.*(1+nu)); % plane strain
            obj.la0 = 2.*obj.mu0.*obj.la0./(obj.la0+2.*obj.mu0); % plane stress 

        end

        function computeStressAndStrain(obj)
            [ux,uy]=pdegrad(obj.p,obj.t,obj.U); % solution gradient
            

            % Seba
            % e=[ux(1,:);(ux(2,:)+uy(1,:))/2;uy(2,:)]; % strain
            % id=[1 0 1]';
            % s=obj.la0*id*(e(1,:)+e(3,:))+2*obj.mu0*e; % stress

            % Swan
            obj.strain = [ux(1,:);(ux(2,:)+uy(1,:))/2;uy(2,:)];
            e = obj.strain;
            obj.tensor = [obj.la0+2*obj.mu0 0 obj.la0; 0 2*obj.mu0 0; obj.la0 0 obj.la0+2*obj.mu0];
            sSwan = obj.tensor*e;
            
            % effective stress
            tgamma3 = [obj.tgamma;obj.tgamma;obj.tgamma]; 
            obj.stress = sSwan.*tgamma3;
        end

        function computeTopologicalDerivative(obj)
            nmat = obj.nMat;
            obj.TD = cell(nmat,nmat);
            E = obj.young;
            a = obj.alpha;
            b = obj.beta;
            s = obj.stress;
            C = obj.tensor;
            %tgamma3 = [obj.tgamma;obj.tgamma;obj.tgamma];
            %e = obj.strain.*tgamma3;

            for i = 1:obj.nMat
                for j = 1:obj.nMat
                    g  = E(j)/E(i);
                    % Seba
                    % coef1 = 0.5*((1-g)./(1+a*g))./obj.tE;
                    % coef2 = coef1.*((g.*(a-2*b)-1)./(1+b*g));
                    % obj.TD{i,j} = coef1.*(4*(s(1,:).*s(1,:)+2*s(2,:).*s(2,:)+s(3,:).*s(3,:))) ...
                    % + coef2.*((s(1,:)+s(3,:)).*(s(1,:)+s(3,:))) ;
                    
                    % Swan
                    c1 = 0.5*((1-g)./(1+a*g))./obj.tE;
                    c2 = c1.*((g.*(a-2*b)-1)./(1+b*g));
                    
                    for k=1:size(s,2)
                        coefMatrix = [4*c1(k)+c2(k) 0 c2(k); 0 8*c1(k) 0; c2(k) 0 4*c1(k)+c2(k)]; 
                        obj.dC(:,:,k,i,j) = C*coefMatrix*C; %(:,:,k,i,j)%
                        %derTop(k) = e(:,k)'*obj.dC(:,:,k,i,j)*e(:,k);
                    end
                    
                    %obj.TD{i,j} = derTop;   
                end      
            end
        end
        

        % function computeVolumeConstraintinDT(obj)
        %     coef = obj.volume(1:end-1) ./ obj.voltarget;
        %     energy = obj.energy0;
        %     pen = obj.penalty;
        %     maxVol = obj.max_vol;
        %     volt = obj.voltarget;
        %     augmentedLagr = obj.auglag;
        %     nmat = obj.nMat;
        % 
        %     for i = 1:obj.nMat
        %         for j = 1:obj.nMat
        %             if obj.penalization == 1
        %                 if i==j
        %                     obj.TD{i,j} = 0;
        %                 elseif j == nmat
        %                     obj.TD{i,j} = obj.TD{i,j}/energy - pen(i)/maxVol;
        %                 elseif i == nmat
        %                     obj.TD{i,j} = obj.TD{i,j}/energy + pen(j)/maxVol;
        %                 else
        %                     obj.TD{i,j} = obj.TD{i,j}/energy + (pen(j) - pen(i))/maxVol;
        %                 end
        %             elseif obj.penalization == 3
        %                 if i==j
        %                     obj.TD{i,j} = 0;
        %                 elseif j == nmat
        %                     obj.TD{i,j} = obj.TD{i,j}/energy - ( (pen(i) + augmentedLagr(i)*(coef(i)-1))/volt(i) );
        %                 elseif i == nmat
        %                     obj.TD{i,j} = obj.TD{i,j}/energy + ( (pen(j) + augmentedLagr(j)*(coef(j)-1))/volt(j) );
        %                 else
        %                     obj.TD{i,j} = obj.TD{i,j}/energy + ( (pen(j) + augmentedLagr(j)*(coef(j)-1))/volt(j) )...
        %                                               - ( (pen(i) + augmentedLagr(i)*(coef(i)-1))/volt(i) );
        %                 end
        %             end
        %         end
        %     end
        % 
        % end
        % 
        % function smoothTopologicalDerivative(obj)
        % 
        %     s.psi = obj.psi;
        %     s.p = obj.p;
        %     s.t = obj.t; 
        %     s.designVariable = obj.designVariable;
        %     s.m = obj.mesh;
        % 
        %     charfun = CharacteristicFunctionComputer(s); 
        %     [~,tfi] = charfun.computeFiandTfi();
        %     [tXi2,~] = integ_exact(obj.t,obj.p,obj.psi(:,2)); chi2 = (1 - tXi2); %- Mixed formulation method
        %     [tXi3,~] = integ_exact(obj.t,obj.p,obj.psi(:,3)); chi3 = (1 - tXi3); %- Mixed formulation method
        %     %     fi = (pdeintrp(p,t,fi)).'; % interpolation at gauss point - P1 projection method
        %     %     chi2 = (pdeintrp(p,t,(psi(:,2)<0))).'; %- P1 projection method
        % 
        %     DT = obj.TD;
        % 
        %     obj.dt = [];
        %     obj.dt(1,:) = - tfi(1,:).*DT{1,end} - tfi(2,:).*DT{2,end} - tfi(3,:).*DT{3,end} ...
        %         + tfi(4,:).*( (1-chi2).*DT{end,1} + (1-chi3).*chi2.*DT{end,2} + chi2.*chi3.*DT{end,3} );
        %     obj.dt(2,:) = - tfi(2,:).*DT{2,1} - tfi(3,:).*DT{3,1} + tfi(1,:).*( (1-chi3).*DT{1,2} + chi3.*DT{1,3} );
        %     obj.dt(3,:) = tfi(2,:).*DT{2,3} - tfi(3,:).*DT{3,2};
        % 
        %     obj.dt = pdeprtni(obj.p,obj.t,obj.dt);
        % end

        
    end




    
        %% material properties
        
        
        %% nominal stress
        
end