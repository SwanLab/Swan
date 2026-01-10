function [NODESV_np1,DATA] = ReactionsFE(BstW,GAUSSV_np1,FORCES,OPERfe,istep,NODESV_np1,wST,Bst,DATA)

if nargin == 0
    load('tmp1.mat')
end

% External forces
% ---------------
F = zeros(size(FORCES{1}.VALUE_ORIGINAL)) ;
for iforces = 1:length(FORCES) % Loop over set of forces (tractions, body ....)
    TIME_FACTOR = FORCES{iforces}.TIME_FACTOR(istep);
    F = F + FORCES{iforces}.VALUE_ORIGINAL*TIME_FACTOR;
end

% Internal forces
if ~isempty(BstW)
Fint = BstW'*GAUSSV_np1.stressST ;
else
    Fint = Bst'*(GAUSSV_np1.stressST.*wST) ;
end

Reactions = zeros(size(F)) ;

%React(DOFr) = K(DOFr,:)*d -F(DOFr) ;
% Slave DOFs
Reactions(OPERfe.DOFs) = Fint(OPERfe.DOFs) - F(OPERfe.DOFs) ;

%  if ~isempty(Gb)
% %         React(DOFm) = - Gb'*React(DOFr) ;
% %     end

if ~isempty(OPERfe.Gbound)
    Reactions(OPERfe.DOFm) = - OPERfe.Gbound'*Reactions(OPERfe.DOFs) ;
end
NODESV_np1.Reactions =  Reactions ; 

DATA = DefaultField(DATA,'ReactionResultantDOFS',[]  ) ; 

if ~isempty(DATA.ReactionResultantDOFS)
    for iii=1:length(DATA.ReactionResultantDOFS)
        
        if  ~isempty(DATA.OPERATORreactions{iii})
            DATA.ReactionResultants(iii,istep)  = DATA.OPERATORreactions{iii}'*(Reactions(DATA.ReactionResultantDOFS{iii}));
        else
            % There might have some problems here... Check versions before
            % 9 M ay to amend it 
            DATA.ReactionResultants(iii,istep)  =   sum(Reactions(DATA.ReactionResultantDOFS{iii}));
        end
    end
end


    


