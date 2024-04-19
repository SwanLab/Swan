matType{1} = load('CircleMicroDamageArea.mat');
matType{2} = load('SquareMicroDamageArea.mat');
matType{3} = load('IsoMicroDamage.mat');

filterTimes1 = 10;
filterTimes2 = 1000;
mat = matType{1}; % 1.Circle, 2.Square, 3.Iso

derivMat = computeGradient(mat,filterTimes1);
deriv2mat = computeGradient(derivMat,filterTimes2);
[funMat,dfunMat,ddfunMat] = computeFunctionsAndDerivatives(mat);

figure()
t = tiledlayout(3,3);
for i=1:3
    for j=1:3
        nexttile
        hold on
        %%%%%%%%%%%%% PLOT DERIVATIVES OF 1 TYPE %%%%%%%%%%%%%%%
        % plot(mat.phi,squeeze(mat.mat(i,j,:)),'Color','#0072BD');
        % plot(derivMat.phi,squeeze(derivMat.mat(i,j,:)),'Color','#D95319');
        % plot(deriv2mat.phi,squeeze(deriv2mat.mat(i,j,:)),'Color','#EDB120');
        % ylabel(['C',num2str(i),num2str(j)]);
        % xlabel("$\phi$",'Interpreter','latex');

        %%%%%%%%%%%%% PLOT ALL TYPES %%%%%%%%%%%%%%%
        % plot(matCircle.phi,squeeze(matCircle.mat(i,j,:)),'Color','#0072BD');
        % plot(matSquare.phi,squeeze(matSquare.mat(i,j,:)),'Color','#D95319');
        % plot(matIso.phi,squeeze(matIso.mat(i,j,:)),'Color','#EDB120');
        % ylabel(['C',num2str(i),num2str(j)]);
        % xlabel("$\phi$",'Interpreter','latex');

        %%%%%%%%%%%%% PLOT DERIVATIVES OF 1 TYPE (POLY VS. PROJECTION) %%%%%%%%%%%%%%%
        % fplot(funMat{i,j},[0 1],'Color','#0072BD');
        % fplot(dfunMat{i,j},[0 1],'Color','#D95319');
        % fplot(ddfunMat{i,j},[0 1],'Color','#EDB120');
        % plot(mat.phi,squeeze(mat.mat(i,j,:)),'.','Color','#0072BD');
        % plot(derivMat.phi,squeeze(derivMat.mat(i,j,:)),'.','Color','#D95319');
        % plot(deriv2mat.phi,squeeze(deriv2mat.mat(i,j,:)),'.','Color','#EDB120');
        % ylabel(['C',num2str(i),num2str(j)]);
        % xlabel("$\phi$",'Interpreter','latex');

        %%%%%%%%%%%%% DERIVATIVES COEFFICIENT ALL TYPES %%%%%%%%%%%%%%%
        for k=1:3
            matTypeInfo = matType{k};
            derivMat = computeGradient(matTypeInfo,filterTimes1);
            deriv2mat = computeGradient(derivMat,filterTimes2);
            coeff = derivMat.mat./deriv2mat.mat;
            plot(matTypeInfo.phi,squeeze(coeff(i,j,:)));
            ylabel('dC/d2C');
            xlabel("$\phi$",'Interpreter','latex');
        end
    end
end
% title(t,['Circle homogenized constitutive tensor (filter ',num2str(filterTimes1),'/',num2str(filterTimes2),')'])
title(t,['Derivatives coefficient (filter ',num2str(filterTimes1),'/',num2str(filterTimes2),')'])
% leg = legend('C','dC','d2C');
leg = legend('Circle','Square','Isotropic');
leg.Layout.Tile = 'east';

function [fun,dfun,ddfun] = computeFunctionsAndDerivatives(cParams)
    x = reshape(cParams.phi,length(cParams.phi),[]);
    y = cParams.mat;
    
    fun   = cell(3,3);
    dfun  = cell(3,3);
    ddfun = cell(3,3);
    for i=1:3
        for j=1:3
            f = fit(x,squeeze(y(i,j,:)),'poly5');
            fun{i,j} = poly2sym(coeffvalues(f));
            dfun{i,j} = diff(fun{i,j});
            ddfun{i,j} = diff(dfun{i,j});
        end
    end
end

function derivMat = computeGradient(cParams,filterTimes)
    x = cParams.phi;
    y = cParams.mat;
    for i=1:3
        for j=1:3
            s.coord = reshape(x,length(x),[]);
            s.connec = (1:1:length(x)-1)' + [0 1];
            mesh1D = Mesh.create(s);
    
            matFun = LagrangianFunction.create(mesh1D,1,'P1');
            matFun.fValues = squeeze(y(i,j,:));
    
            derivMatDom = Grad(matFun);
            derivMatFun = derivMatDom.project('P1',mesh1D);
            if filterTimes ~= 0
                for k=1:filterTimes
                    ss.filterType = 'LUMP';
                    ss.mesh       = mesh1D;
                    ss.trial      = LagrangianFunction.create(mesh1D,1,'P1');
                    filter        = Filter.create(ss);
                    derivMatFun   = filter.compute(derivMatFun,'QUADRATIC');
                end  
            end


            derivMat.mat(i,j,:) = derivMatFun.fValues;
            
        end
    end
    derivMat.phi = cParams.phi;
end