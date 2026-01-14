function  VAR =  ResidualDynamicPart(VAR,OPERFE,TIMEPARAM,DeltaT)

if nargin == 0
    load('tmp2.mat')
end



if TIMEPARAM.TYPE==1
    % Newmark-Bossak
    
    %  \begin{equation}
    %  \label{eq:***1s}
    %  \begin{split}
    %    \Res{}{}(\d) & =   \Par{  \{\FintNEWMARK{}{}(\d)}  + \Fint{}{}(\d)  }
    %      - \Par{ \{\FextNEWMARK{}{n+1}}   +   \Fext{}{n+1}}
    %    \end{split}
    %  \end{equation}
    
    % while the inertial/damping forces by
    %  \begin{equation}
    %  \label{eq:::--:0}
    %  \cblue{ \FintNEWMARK{}{}(\d) \defeq  \cblue{\Par{\dfrac{1-\aBOSSAK}{\beta \Delta t_{n+1}^2} \M_{}  +  \dfrac{\gamma}{\beta \Delta t_{n+1}} \Ddamp_{} } \d} }
    %  \end{equation}
    
    % Contribution displacement-independent damp./iner. external forces
    VAR.RESID = VAR.RESID - VAR.FextENB ;
    
    % Displacement-dependent term
    % ------------------------------
    FintNEWMARK = (1-TIMEPARAM.a)/(TIMEPARAM.BETA*DeltaT^2)*OPERFE.M*VAR.DISP  ;
    if ~isempty(OPERFE.Ddamp)
        FintNEWMARK = FintNEWMARK + TIMEPARAM.GAMMA/(TIMEPARAM.BETA*DeltaT)*OPERFE.Ddamp*VAR.DISP ;
    end
    VAR.RESID = VAR.RESID + FintNEWMARK ; 
    
    
end