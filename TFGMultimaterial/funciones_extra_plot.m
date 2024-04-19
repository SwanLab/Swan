
% %For a comparison between shape funcion and potential energy
%  figure(4); clf;hold on, plot(gsf,'k-','LineWidth',1.2); %title('Shape Function');
%  yyaxis left; ylabel('Shape Function','Interpreter','latex','FontSize',18);
%  yyaxis right; plot(gEpot,'-','LineWidth',1.2);   ylabel('Epot','Interpreter','latex','FontSize',18);
%  grid minor, xlim([0 iter]), ylim([min([gsf,gEpot]) max([gsf,gEpot])])
%   ax = gca; ax.TickLabelInterpreter = 'latex'; ax.FontSize = 16;
  
  %For a comparison between shape funcion and potential energy -- V2
%   Color = [0 113 188]/255; %2 mat
%   Color = [216 82 24]/255; %3 mat
  Color = [236 176 31]/255; %4 mat
  figure(10); clf; plot(gsf,'-','Color',Color,'LineWidth',1.2); hold on, %title('Shape Function');
  plot(gEpot,'--','Color',Color,'LineWidth',1.2);
%   plot(gL,'--','Color',Color,'LineWidth',1.2);
  grid minor, xlim([0 iter]),
  ax = gca; ax.TickLabelInterpreter = 'latex'; ax.FontSize = 16;
  legend({'$\mathcal{J}(\chi)$','$\mathcal{L}(\chi)$'},'Interpreter','latex','FontSize',14,'Location','southeast')
  xlabel('Iterations','Interpreter','latex','FontSize',18)

% figure(4); clf; semilogx(gsf,'LineWidth',1.2); hold on, 
% ylabel('Shape Function','Interpreter','latex','FontSize',18);
% xlabel('Iter','Interpreter','latex','FontSize',16);
% legend({'1 mat','2 mat','3 mat'},'Location','northeast','Interpreter','latex','FontSize',14),
% ax = gca; ax.TickLabelInterpreter = 'latex'; ax.FontSize = 16;
% grid minor,
 
 figure(5); clf; plot(gth*180/pi,'-','MarkerFaceColor',facecolor,'Color',linecolor,'LineWidth',1.2); title('Theta Angle');
 grid minor , xlim([0 iter]),ylim([min(gth*180/pi)-2 max(gth*180/pi)+10])
 drawnow
 
 figure(6); clf;
 for index = 1:size(gvol,1)
     plot(gvol(index,:),'-','LineWidth',1.2), hold on
     leyenda{index} = sprintf('$V_{f%d}$',index);
 end
 legend(leyenda,'Interpreter','latex'),
 ax = gca; ax.TickLabelInterpreter = 'latex'; ax.FontSize = 16;
 ylabel('Volume Fraction','Interpreter','latex','FontSize',18);
 xlabel('Iterations','Interpreter','latex','FontSize',18)
 grid minor, 
 xlim([0 iter]), ylim([-5 105]),
 leyenda{end+1} = 'target';
 hold on, plot([0 iter],params.volfrac(1)*100*[1 1],'k--'),legend(leyenda)


 figure(20),clf,
 maxIter = max([length(gth_1),length(gth_2),length(gth_3)]);
 semilogx(gth_1*180/pi), hold on,semilogx(gth_2*180/pi), hold on,semilogx(gth_3*180/pi),
 semilogx([1 maxIter],[1 1],'k--'),  hold on,xLIM = get(gca,'XLim');
 ylim([-1 inf]),xlim([0 maxIter])
 legend({'1 mat','2 mat','3 mat','$\theta_{stop}$'},'Location','northeast','Interpreter','latex','FontSize',14), grid minor,
 ylabel('$\theta$','Interpreter','latex','FontSize',18)
 xlabel('Iterations','Interpreter','latex','FontSize',18)
 set(gca,'TickLabelInterpreter','latex','FontSize',16)
 axes('Position',[0.495 0.31 0.4 0.35])
 box on
 xzoom = 15; 
 iter_1 = 1:length(gth_1); index_1=(iter_1>=xzoom & iter_1<=maxIter);
 iter_2 = 1:length(gth_2); index_2=(iter_2>=xzoom & iter_2<=maxIter);
 iter_3 = 1:length(gth_3); index_3=(iter_3>=xzoom & iter_3<=maxIter);
 semilogx([xzoom maxIter],[1 1],'k--'),  hold on,
 semilogx(iter_1(index_1),gth_1(index_1)*180/pi), hold on, semilogx(iter_2(index_2),gth_2(index_2)*180/pi), hold on, semilogx(iter_3(index_3),gth_3(index_3)*180/pi), 
 hold on, xlim([xzoom maxIter+5])
 grid minor, set(gca,'TickLabelInterpreter','latex','FontSize',12)

 %For ploting images matrices and mirrors

 %Get rid of background
% HH = figure('Menu','none','ToolBar','none'); 
HH = figure('ToolBar','none'); 
clf;
ah = axes('Units','Normalize','Position',[0 0 1 1]);
        multimat_plot( p,t,fi );
        axis(ah, 'square')
HH.Position(3) = HH.Position(4);
% print (HH, '-deps',  strcat(filename,'_topology'));
nuevoframe = getframe(HH);

%Upload image
% imagen = imread('symmetry_cog_3m_topology.png'); 
imagen = nuevoframe.cdata;

%For making matrices of images
%Opcion 1
% imagenArray = repmat(imagen, [ 1 1 1 4 ]);
% montage(imagenArray);
%Opcion 2
imagenArray = repmat(imagen, 50);
imshow(imagenArray);

%For making mirrors
imagen2 = flip(imagen,2);
imagen12 = [imagen2 imagen];
imagen3 = flip(imagen12,1);
imagen4 = [imagen12; imagen3];
% imshow(imagen);
% imshow(imagen2);
% imshow(imagen3);

fig = figure('ToolBar','none'); 
clf;
ax = axes('Units','Normalize','Position',[0 0 1 1]);
        imshow(imagen4);
        axis(ax, 'square')
fig.Position(3) = fig.Position(4);
 
 