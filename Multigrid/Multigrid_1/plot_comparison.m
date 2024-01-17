% %Post-processing
subplot(2,2,1)
surf(x_mesh,y_mesh,trap_sol_exp)
xlabel('X');
ylabel('Y');
title('FDM 5 point stencil Implicit')
axis equal;
colorbar;
caxis([-0.5, 0.5]);
colorbar('off')

subplot(2,2,2)
surf(x_mesh,y_mesh,imp_sol)
xlabel('X');
ylabel('Y');
title('FEM explicit')
axis equal;
colorbar('Ticks', linspace(-0.5,.5,11));
caxis([-0.5, 0.5]);

error_mesh = trap_sol_exp - imp_sol;
subplot(2,2,[3,4]);
surf(x_mesh,y_mesh,error_mesh)
xlabel('X');
ylabel('Y');
colorbar;
title('FDM - FEM')
axis equal;
str = sprintf('h = %4.3f (%d elements), dt = %4.3f ', h, nelem, dt);
suptitle(strcat(str, ': Comparison of explicit FEM and implicit FDM'));
    
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 26 20])
set(gcf,'PaperOrientation','landscape');
print(gcf, '-dpdf', 'FDM_FEM_exp.pdf')