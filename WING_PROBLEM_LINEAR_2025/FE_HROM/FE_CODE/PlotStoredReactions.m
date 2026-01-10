function ReactionsFINAL = PlotStoredReactions(DATA_nameWORKSPACE)

 ReactionsFINAL = [] ; 

load(DATA_nameWORKSPACE,'ReactionsPlot','ReactionsFINAL') ;

figure(ReactionsPlot.ifigure) ;
hold on
grid on
xlabel('Time factor') ;
ylabel('Force ') ;
grid on
title(['Reaction forces '])
hplot = []   ;
for iii = 1:length(ReactionsPlot.xvar)
%     if ~isempty(DATA)
%         x = DATA.FACTOR_TIME_DISPLACEMENTS ; 
%         if length(x) ~= length(ReactionsPlot.yvar{iii})  ;% To amend an error detected in Jan-23-2020
%             x = ReactionsPlot.xvar{iii} ;
%         end
%     end
    hplot(end+1) = plot(ReactionsPlot.xvar{iii},ReactionsPlot.yvar{iii},ReactionsPlot.COLORS_STORE{iii} ) ;
end

legend(hplot,ReactionsPlot.LEGG) ;