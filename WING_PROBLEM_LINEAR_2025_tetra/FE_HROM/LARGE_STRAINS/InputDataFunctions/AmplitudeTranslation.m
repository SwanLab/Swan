function AMPLITUDE = AmplitudeTranslation(TRANSLATION,INTERVAL,TIMELOC,rnodLOC)

FactorSteps =  TRANSLATION.TIMEFUN(TIMELOC) ;
FactorSteps = FactorSteps.*(TIMELOC >= INTERVAL(1)) ;
FactorSteps = FactorSteps.*(TIMELOC <= INTERVAL(2)) ;
AMPLITUDE= FactorSteps*TRANSLATION.AMPLITUDE(:) ;
AMPLITUDE = repmat(AMPLITUDE,length(rnodLOC),1) ;

end

