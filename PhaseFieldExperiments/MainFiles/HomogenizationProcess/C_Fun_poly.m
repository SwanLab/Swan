%%%%%%%%%%%%%%%%%%%%%%%%%%% ALL TYPE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matType{1} = load('CircleMicroDamageArea.mat');
matType{2} = load('CircleMicroDamagePerimeter.mat');
matType{3} = load('SquareMicroDamageArea.mat');
matType{4} = load('SquareMicroDamagePerimeter.mat');
matType{5} = load('IsoMicroDamage.mat');

filterTimes1 = 0;
filterTimes2 = 0;
 for i=1:length(matType)
%     derivMatType{i} = computeGradient(matType{i},filterTimes1);
%     deriv2MatType{i} = computeGradient(derivMatType{i},filterTimes2);
     [funMat(:,:,i),dfunMat(:,:,i),ddfunMat(:,:,i)] = computeFunctionsAndDerivatives(matType{i});
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%% ONE TYPE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

type = 'Circle (Perimeter)';
switch type
    case 'Circle (Area)'
        mat = matType{1};
        derivMat = derivMatType{1};
        deriv2Mat = deriv2MatType{1};
    case 'Circle (Perimeter)'
        mat = matType{2};
%        derivMat = derivMatType{2};
%        deriv2Mat = deriv2MatType{2};
    case 'Square (Area)'
        mat = matType{3};
        derivMat = derivMatType{3};
        deriv2Mat = deriv2MatType{3};
    case 'Square (Perimeter)'
        mat = matType{4};
        derivMat = derivMatType{4};
        deriv2Mat = deriv2MatType{4};
    case 'Isotropic'
        mat = matType{5};
        derivMat = derivMatType{5};
        deriv2Mat = deriv2MatType{5};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure()
t = tiledlayout(3,3);
for i=1:3
    for j=1:3
        nexttile
        hold on
        %%%%%%%%%%%%% PLOT DERIVATIVES OF 1 TYPE %%%%%%%%%%%%%%%
        % plot(mat.phi,squeeze(mat.mat(i,j,:)),'Color','#0072BD');
        % plot(derivMat.phi,squeeze(derivMat.mat(i,j,:)),'Color','#D95319');
        % plot(deriv2Mat.phi,squeeze(deriv2Mat.mat(i,j,:)),'Color','#EDB120');
        % ylabel(['C',num2str(i),num2str(j)]);
        % xlabel("$\phi$",'Interpreter','latex');
        % title(t,[type,' homogenized constitutive tensor (filter ',num2str(filterTimes1),'/',num2str(filterTimes2),')'])
        % legType = 1;

        %%%%%%%%%%%%% PLOT TENSOR ALL TYPES %%%%%%%%%%%%%%%
        % plot((matType{1}.phi),squeeze(matType{1}.mat(i,j,:)),'Color','#0072BD');
        % plot((matType{2}.phi),squeeze(matType{2}.mat(i,j,:)),'--','Color','#0072BD');
        % plot((matType{3}.phi),squeeze(matType{3}.mat(i,j,:)),'Color','#D95319');
        % plot((matType{4}.phi),squeeze(matType{4}.mat(i,j,:)),'--','Color','#D95319');
        % plot(sqrt((matType{5}.phi)),squeeze(matType{5}.mat(i,j,:)),'Color','#EDB120');
        % plot((matType{5}.phi),squeeze(matType{5}.mat(i,j,:)),'k');
        % plot((matType{5}.phi).^2,squeeze(matType{5}.mat(i,j,:)),'--','Color','#EDB120');
        % ylabel(['C',num2str(i),num2str(j)]);
        % xlabel("$\phi$",'Interpreter','latex');
        % title(t,'Homogenized constitutive tensors')
        % legType = 2;

        %%%%%%%%%%%%% PLOT DERIVATIVE TENSOR ALL TYPES %%%%%%%%%%%%%%%
        % plot(derivMatType{1}.phi,squeeze(derivMatType{1}.mat(i,j,:)),'.','Color','#0072BD');
        % plot(derivMatType{2}.phi,squeeze(derivMatType{2}.mat(i,j,:)),'Color','#0072BD');
        % plot(derivMatType{3}.phi,squeeze(derivMatType{3}.mat(i,j,:)),'.','Color','#D95319');
        % plot(derivMatType{4}.phi,squeeze(derivMatType{4}.mat(i,j,:)),'Color','#D95319');
        % plot(derivMatType{5}.phi,squeeze(derivMatType{5}.mat(i,j,:)),'Color','#EDB120');
        % 
        % ylabel(['C',num2str(i),num2str(j)]);
        % xlabel("$\phi$",'Interpreter','latex');
        % title(t,'Homogenized constitutive tensors derivatives')
        % legType = 2;

        %%%%%%%%%%%%% PLOT SECOND DERIVATIVE TENSOR ALL TYPES %%%%%%%%%%%%%%%
        % plot(deriv2MatType{1}.phi,squeeze(deriv2MatType{1}.mat(i,j,:)),'.','Color','#0072BD');
        % plot(deriv2MatType{2}.phi,squeeze(deriv2MatType{2}.mat(i,j,:)),'Color','#0072BD');
        % plot(deriv2MatType{3}.phi,squeeze(deriv2MatType{3}.mat(i,j,:)),'.','Color','#D95319');
        % plot(deriv2MatType{4}.phi,squeeze(deriv2MatType{4}.mat(i,j,:)),'Color','#D95319');
        % plot(deriv2MatType{5}.phi,squeeze(deriv2MatType{5}.mat(i,j,:)),'Color','#EDB120');
        % 
        % ylabel(['C',num2str(i),num2str(j)]);
        % xlabel("$\phi$",'Interpreter','latex');
        % title(t,'Homogenized constitutive tensors derivatives')
        % legType = 2;

        %%%%%%%%%%%%% PLOT DERIVATIVES OF 1 TYPE (POLY VS. PROJECTION) %%%%%%%%%%%%%%%
        % fplot(funMat{i,j},[0 1],'Color','#0072BD');
        % fplot(dfunMat{i,j},[0 1],'Color','#D95319');
        % fplot(ddfunMat{i,j},[0 1],'Color','#EDB120');
        % plot(mat.phi,squeeze(mat.mat(i,j,:)),'.','Color','#0072BD');
        % %plot(derivMat.phi,squeeze(derivMat.mat(i,j,:)),'.','Color','#D95319');
        % %plot(deriv2Mat.phi,squeeze(deriv2Mat.mat(i,j,:)),'.','Color','#EDB120');
        % ylabel(['C',num2str(i),num2str(j)]);
        % xlabel("$\phi$",'Interpreter','latex');
        % title(t,[type,' homogenized constitutive tensor (filter ',num2str(filterTimes1),'/',num2str(filterTimes2),')'])
        % legType = 1;

        %%%%%%%%%%%%%%%%%%% PLOT ANALYTICAL ALL TYPES %%%%%%%%%%%%%%%%%%
        % fplot(funMat(i,j,1),[0 1],'--','Color','#0072BD');
        % fplot(funMat(i,j,2),[0 1],'Color','#0072BD');
        % fplot(funMat(i,j,3),[0 1],'--','Color','#D95319');
        % fplot(funMat(i,j,4),[0 1],'Color','#D95319');
        % fplot(funMat(i,j,5),[0 1],'Color','#EDB120');
        % ylabel(['C',num2str(i),num2str(j)]);
        % xlabel("$\phi$",'Interpreter','latex');
        % title(t,'ANALYTICAL homogenized constitutive tensor')
        % legType = 2;

        %%%%%%%%%%%%%%%%%%%% PLOT DERIVATIVES ALL TYPES %%%%%%%%%%%%%%%%%
        % fplot(dfunMat(i,j,1),[0 1],'Color','#0072BD');
        % fplot(dfunMat(i,j,2),[0 1],'--','Color','#0072BD');
        % fplot(dfunMat(i,j,3),[0 1],'Color','#D95319');
        % fplot(dfunMat(i,j,4),[0 1],'--','Color','#D95319');
        % fplot(dfunMat(i,j,5),[0 1],'Color','#EDB120');
        % ylabel(['C',num2str(i),num2str(j)]);
        % xlabel("$\phi$",'Interpreter','latex');
        % title(t,'ANALYTICAL homogenized constitutive tensor DERIVATIVE')
        % legType = 2;

        %%%%%%%%%%%%%%%%%%% PLOT SECOND DERIVATIVES ALL TYPES %%%%%%%%%%%%%
        % fplot(ddfunMat(i,j,1),[0 1],'Color','#0072BD');
        % fplot(ddfunMat(i,j,2),[0 1],'--','Color','#0072BD');
        % fplot(ddfunMat(i,j,3),[0 1],'Color','#D95319');
        % fplot(ddfunMat(i,j,4),[0 1],'--','Color','#D95319');
        % fplot(ddfunMat(i,j,5),[0 1],'Color','#EDB120');
        % ylabel(['C',num2str(i),num2str(j)]);
        % xlabel("$\phi$",'Interpreter','latex');
        % title(t,'ANALYTICAL homogenized constitutive tensor SECOND DERIVATIVE')
        % legType = 2;

        %%%%%%%%%%%%% DERIVATIVES COEFFICIENT ALL TYPES %%%%%%%%%%%%%%%
        % for k=1:length(matType)
        %     matTypeInfo = matType{k};
        %     derivMat = computeGradient(matTypeInfo,filterTimes1);
        %     deriv2Mat = computeGradient(derivMat,filterTimes2);
        %     coeff = derivMat.mat./deriv2Mat.mat;
        %     plot(matTypeInfo.phi,squeeze(coeff(i,j,:)));
        %     ylabel('dC/d2C');
        %     xlabel("$\phi$",'Interpreter','latex');
        % end
        % title(t,['Derivatives coefficient (filter ',num2str(filterTimes1),'/',num2str(filterTimes2),')'])
        % legType = 2;

    end
end

switch legType
    case 1
        leg = legend('C','dC','d2C');
    case 2
        leg =legend('Circle (Area)','Circle (Perimeter)','Square (Area)','Square (Perimeter)','AT0.5','AT1','AT2');
    case 3
        leg = legend('Circle (Perimeter)','Square (Perimeter)','Isotropic');
end
leg.Layout.Tile = 'east';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [fun,dfun,ddfun] = computeFunctionsAndDerivatives(cParams)
    x = reshape(cParams.phi,length(cParams.phi),[]);
    y = cParams.mat;
    
    fun   = cell(3,3);
    dfun  = cell(3,3);
    ddfun = cell(3,3);
    for i=1:3
        for j=1:3
            f = fit(x,squeeze(y(i,j,:)),'poly9');
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