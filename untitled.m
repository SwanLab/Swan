matCircle = load('CircleMicroDamagePerimeter.mat');
matSquare = load('SquareMicroDamagePerimeter.mat');
matIso    = load('IsoMicroDamage.mat');

filterTimes1 = 0;
filterTimes2 = 0;
mat = matIso;
derivMat = computeGradient(mat,filterTimes1);
deriv2mat = computeGradient(derivMat,filterTimes2);
close all
t = tiledlayout(3,3);
for i=1:3
    for j=1:3
        nexttile
        hold on

        % plot(mat.phi,squeeze(mat.mat(i,j,:)),'Color','#0072BD');
        % plot(derivMat.phi,squeeze(derivMat.mat(i,j,:)),'Color','#D95319');
        % plot(deriv2mat.phi,squeeze(deriv2mat.mat(i,j,:)),'Color','#EDB120');

        plot(matCircle.phi,squeeze(matCircle.mat(i,j,:)),'Color','#0072BD');
        plot(matSquare.phi,squeeze(matSquare.mat(i,j,:)),'Color','#D95319');
        plot(matIso.phi,squeeze(matIso.mat(i,j,:)),'Color','#EDB120');

        ylabel(['C',num2str(i),num2str(j)]);
        xlabel("$\phi$",'Interpreter','latex');
    end
end

% title(t,['Isotropic constitutive tensor (filter ',num2str(filterTimes1),'/',num2str(filterTimes2),')'])
% leg = legend('C','dC','d2C');

title(t,'Material constitutive tensor (damage = Perimeter)');
leg = legend('Circle','Square','Isotropic');
leg.Layout.Tile = 'east';


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
                    derivMatFun     = filter.compute(derivMatFun,'QUADRATIC');
                end  
            end


            derivMat.mat(i,j,:) = derivMatFun.fValues;
            
        end
    end
    derivMat.phi = cParams.phi;
end



% s0.coord = (x(2:end)+x(1:end-1))/2;
% s0.connec = (1:1:length(x)-2)' + [0 1];
% meshP0 = Mesh.create(s0);
% 
% P1 = LagrangianFunction.create(mesh1D,1,'P1');
% P1.fValues = y;
% figure()
% P1.plot
% title('P1')
% 
% GP1Dom = Grad(P1);
% GP1 = GP1Dom.project('P1',mesh1D);
% 
% GP1Deriv = LagrangianFunction.create(meshP0,1,'P1');
% GP1Deriv.fValues = squeeze(GP1Dom.evaluate(0));
% for i=1:10
%     ss.filterType = 'LUMP';
%     ss.mesh       = meshP0;
%     ss.trial      = LagrangianFunction.create(meshP0,1,'P1');
%     filter        = Filter.create(ss);
%     GP1Deriv      = filter.compute(GP1Deriv,'QUADRATIC');
% end   
% figure()
% GP1Deriv.plot
% title('GP1')
% 
% 
% 
% GGP1Dom = Grad(GP1);
% GGP1 = GGP1Dom.project('P1',mesh1D);
% figure()
% GGP1.plot
% title('GGP1')
% for i=1:10
%     ss.filterType = 'LUMP';
%     ss.mesh       = mesh1D;
%     ss.trial      = LagrangianFunction.create(mesh1D,1,'P1');
%     filter        = Filter.create(ss);
%     GGP1       = filter.compute(GGP1,'QUADRATIC');
% end
% figure()
% title('Figure 11')
% GGP1.plot()
% 
% 
% % ss.filterType   = 'PDE';
% % ss.mesh         = mesh1D;
% % ss.boundaryType = 'Robin';
% % ss.metric       = 'Isotropy';
% % filter          = Filter.create(ss);
% % 
% % epsilon    = 4*mesh1D.computeMeanCellSize();
% % filter.updateEpsilon(epsilon);
% % fP1PDE = filter.compute(GGP1,'QUADRATIC');
% % figure()
% % title('Figure 8')
% % fP1PDE.plot()
% 
% P2 = P1.project('P2');
% figure()
% P2.plot
% title('P2')
% 
% GP2Dom = Grad(P2);
% GP2 = GP1Dom.project('P1',mesh1D);
% figure()
% GP2.plot
% title('GP2')
% 
% GGP2Dom = Grad(GP2);
% GGP2 = GGP2Dom.project('P2',mesh1D);
% figure()
% GGP2.plot
% title('GGP2')
% 
% y = P2.evaluateHessian(0);
% figure()
% plot(s0.coord,squeeze(y))