function plotErrors(DispError,StressError,numberElements)
    RelativeTipError = DispError.relative;
    RelativeStressError = StressError.relative;

    figure
    loglog(numberElements,squeeze(RelativeTipError),'-o')
    grid on
    title('Relative error vs. elements')
    xlabel('Number of elements')
    ylabel('Relative error')
    ytickformat('%.2f')
    legend('Tip displacement (x=1,y=0.25)')

    figure
    loglog(numberElements,squeeze(RelativeStressError(1,1,:)),'-o')
    hold on
    loglog(numberElements,squeeze(RelativeStressError(1,2,:)),'-o')
    loglog(numberElements,squeeze(RelativeStressError(2,1,:)),'-o')
    loglog(numberElements,squeeze(RelativeStressError(2,2,:)),'-o')
    grid on
    title('Relative error vs. elements')
    xlabel('Number of elements')
    ylabel('Relative error')
    ytickformat('%.3f')
    legend('Normal stress (x=0,y=0.05)', ...
            'Normal stress (x=0,y=0.1)', ...
            'Shear stress (x=0,y=0.05)', ...
            'Shear stress (x=0,y=0.1)')
end