function DATA  = ComputeReactionsONLINE(DATA,OPERfe)
if nargin == 0
    load('tmp.mat')
end

DATA = DefaultField(DATA,'IS_MULTISCALE_ROM_ELEMENT',0) ; 

if length(OPERfe.ndim) == 1 && DATA.IS_MULTISCALE_ROM_ELEMENT ==0  % Standard element
    DATA  = ReactionPlotONLINE_standard(DATA,OPERfe) ;
elseif length(OPERfe.ndim) > 1
    % Multiscale element
    DATA  = ReactionPlotONLINE_multiscale(DATA,OPERfe) ;
  
elseif DATA.IS_MULTISCALE_ROM_ELEMENT ==1
    DATA  = ReactionPlotONLINE_multiscale1D(DATA,OPERfe) ;
end
