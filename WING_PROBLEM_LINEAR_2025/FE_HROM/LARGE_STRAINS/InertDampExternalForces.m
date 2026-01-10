function   VAR = InertDampExternalForces(VAR,TIME_INTEGRATION_PARAMETERS,DeltaT,OPERFE)

if nargin == 0
    load('tmp2.mat')
end




if TIME_INTEGRATION_PARAMETERS.TYPE==1
    % Bossak-Newmark integration scheme
    % ---------------------------------
    ALPHA = TIME_INTEGRATION_PARAMETERS.a ;
    BETA = TIME_INTEGRATION_PARAMETERS.BETA ;
    GAMMA = TIME_INTEGRATION_PARAMETERS.GAMMA ; 
    % \begin{equation*}
    % \label{eq::987}
    %  \cblue{\dTIL^{n+1} = \d^n + \Delta t_{n+1} \v^n  + \dfrac{\Delta t_{n+1}^2}{2} (1-2 \beta ) \a^n}
    % \end{equation*}
    VAR.DISPtilde =VAR.DISP +  DeltaT*VAR.VEL + 0.5*DeltaT^2*(1-2*BETA)*VAR.ACEL ;
    %
    %\begin{equation*}
    % {\vTIL^{n+1} = \v^n + (1-\gamma) \Delta t_{n+1} \a^n }
    %\end{equation*}
    VAR.VELtilde = VAR.VEL + (1-GAMMA)*DeltaT*VAR.ACEL ;
    
    % External forces
    
%     \begin{equation*}
%  \label{eq:5670*}
%  \begin{split}
%    \FextNEWMARK{}{\DOFl} &    \M_{\DOFl} (  \dfrac{1}{\beta \Delta t_{n+1}^2} \dTIL^{n+1} - \aBOSSAK \a^n  ) + \Ddamp_{\DOFl}   \Par{ \dfrac{\gamma}{\beta \Delta t_{n+1}} \dTIL^{n+1} -   \vTIL^{n+1}  }
%   \end{split} 
%  \end{equation*}


    VAR.FextENB = OPERFE.M*((1-ALPHA)/BETA/DeltaT^2*VAR.DISPtilde - ALPHA*VAR.ACEL )  ; 
    if ~isempty(OPERFE.Ddamp) 
        VAR.FextENB = VAR.FextENB + OPERFE.Ddamp*(GAMMA/(BETA*DeltaT)*VAR.DISPtilde - VAR.VELtilde)  ; 
    end
    
else
    error('Option not implemented')
end
