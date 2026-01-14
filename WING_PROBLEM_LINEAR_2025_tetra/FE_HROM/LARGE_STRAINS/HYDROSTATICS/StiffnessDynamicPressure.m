function kdyn = StiffnessDynamicPressure(k_m,b,m,DATA,v,vn_pos,n,OPERHYDRO,densW)
%   \begin{equation}
%  \label{eq:4.m}
%  \k_{dyn} =    \densW  v_n  \cblue{\derpar{v_n}{\d^e} }\HEAV{v_n} \HEAV{y}
%  \end{equation}
% B = der(vn)
if nargin  == 0
    load('tmp.mat')
end
% b = norm(m)
% k_m = der(m)/der(d^e)
% Velocity
% Normal velocity
% -----------------------------
% Output \derpar{v_n}{d_k}

% ---------------------------------------------------------------------------------------------------------------
% % \begin{equation}
%   \derpar{v_n}{d_k} =   \cblue{\derpar{n_i}{d_k}} v_i   +  n_i  N_{ij} \cred{\derpar{\dot{d}_j }{d_k}}
%  \end{equation}
% ---------------------------------------------------------------------------------------------------------------
%  B_k =  \ {\derpar{n_i}{d_k}} v_i

%  C_k = n_i  N_{ij} \cred{\derpar{\dot{d}_j }{d_k}}


% -----------------------------------------------------------------------------------------------------------------
%  \item Term $\cblue{\derpar{n_i}{d_k}}$.  Since
%   \begin{equation}
%   \boxed{\cblue{\derpar{n_i}{d_k}}  = \cgreen{\derpar{n_i}{m_a}}  \cred{\derpar{m_a}{d_k}}}
%  \end{equation}
% ------------------------------------------------------------------------------------------------------------------




% TErm  ---  G_{ia} \defeq \cgreen{\derpar{n_i}{m_a}}  =   \dfrac{\delta_{ia}  b^2  -  m_i m_a }{b^3}}
ndim  = DATA.MESH.ndim ;
G = zeros(size(m,1),ndim) ;

ngausTOT= length(m)/ndim ;
igaus = 1:ngausTOT ;
for idim = 1:ndim
    idimGLO = (igaus-1)*ndim + idim ;
    for adim = 1:ndim
        adimGLO = (igaus-1)*ndim + adim ;
        if idim == adim
            G(idimGLO,adim)  =  b.^2 ;
        end
        G(idimGLO,adim)  =  G(idimGLO,adim)  - m(idimGLO).*m(adimGLO) ;
        G(idimGLO,adim) =  G(idimGLO,adim)./(b.^3)  ;
    end
end






% ---- \cred{\derpar{m_a}{d_k}}}  ---->  k_m
% --------------------------------------------
% Therefore
% B_k =  \boxed{\cblue{\derpar{n_i}{d_k}}*v_i  = \cgreen{\derpar{n_i}{m_a}}
% \cred{\derpar{m_a}{d_k}}}*v_i  = G_ia (k_m)_ak*v_i

%   B =  \derpar{n_i}{d_k}*v_i
B = zeros(length(b),size(k_m,2));



VECTmax = 1 ;  % Optimized

if VECTmax == 0
    % Not used, just for trials
    
    for  idim = 1:ndim
        idimGLO = (igaus-1)*ndim + idim ;
        for kdim = 1:size(B,2)
            for adim = 1:ndim
                adimGLO = (igaus-1)*ndim + adim ;
                B(:,kdim)  =  B(:,kdim) + G(idimGLO,adim).*k_m(adimGLO,kdim).*v(idimGLO) ;
            end
        end
    end
    
else
    
    for  idim = 1:ndim
        idimGLO = (igaus-1)*ndim + idim ;
        %    for kdim = 1:size(B,2)
        kdim = 1:size(B,2) ;
        for adim = 1:ndim
            adimGLO = (igaus-1)*ndim + adim ;
            B(:,kdim)  =  B(:,kdim) +  bsxfun(@times,k_m(adimGLO,kdim),G(idimGLO,adim).*v(idimGLO)) ; % G(idimGLO,adim).*k_m(adimGLO,kdim).*v(idimGLO) ;
        end
        %       end
    end
    
    
end



%%% LEt's attack the second term:   %  C_k = n_i  N_{ij}
%%% \cred{\derpar{\dot{d}_j }{d_k}}
% ----------------------------------------------------------------------------------------------

% \derpar{\dot{d}_j }{d_k}

%   \begin{equation}
%    \derpar{\v}{\d}   =   = \dfrac{\gamma  }{ \beta \Delta t} \ident
%   \end{equation}

% Therefore
% %  \begin{equation}
%    n_i  N_{ij} \cred{\derpar{\dot{d}_j }{d_k}}  = n_i N_{ij}  \delta_{jk}   = \dfrac{\gamma  }{ \beta \Delta t}     n_i N_{ik}
%   \end{equation}

if DATA.TIME_INTEGRATION_PARAMETERS.TYPE==1
    % Bossak-Newmark integration scheme
    % ---------------------------------
    
    BETA = DATA.TIME_INTEGRATION_PARAMETERS.BETA ;
    GAMMA = DATA.TIME_INTEGRATION_PARAMETERS.GAMMA ;
    cte = GAMMA/BETA/DATA.DeltaT ;
    
    
    VECTmax = 1 ;  % Optimized

    if VECTmax == 0
        
        
        for idim = 1:ndim
            idimGLO = (igaus-1)*ndim + idim ;
            for kdim = 1:size(B,2)
                B(:,kdim) = B(:,kdim)  + cte*n(idimGLO).*OPERHYDRO.NbST(idimGLO,kdim) ;
            end
        end
        
    else
        for idim = 1:ndim
            idimGLO = (igaus-1)*ndim + idim ;
            kdim = 1:size(B,2)  ;
            B(:,kdim) = B(:,kdim)  + cte*bsxfun(@times,OPERHYDRO.NbST(idimGLO,kdim),n(idimGLO))   ;
            
        end
        
    end
    
else
    error('Option not implemented')
end

% %  \k_{dyn} =    \densW  v_nPOS  \cblue{\derpar{v_n}{\d^e} }

kdyn = densW*bsxfun(@times,B,vn_pos) ;
