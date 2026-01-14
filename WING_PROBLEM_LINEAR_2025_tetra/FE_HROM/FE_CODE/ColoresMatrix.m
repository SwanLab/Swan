function colores = ColoresMatrix(ntrial,nonrep,varargin);
% Matriz de colores
% For avoiding white  --> ntrial+1


if nargin == 1
    nonrep = [] ;
elseif nargin == 0
    nonrep = [] ;
    ntrial = 1 ;
end

% LOCAL INPUTS
ALEATORIO = 'YES';
% END LOCAL INPUTS

%% EXTRACTING INPUTS
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varginOR = varargin ;
FdnamesInputs = {'ALEATORIO'};
AuxDATA=WriteAuxFdNamesNEW(FdnamesInputs,varginOR);
for id = 1:length(AuxDATA);
    eval(AuxDATA{id});
end




yellow = [1 1 0];
ntrialor = ntrial ;

if isempty(nonrep)
    nonrep = yellow; % No yellow
else
    nonrep = [nonrep;1 1 0];



end
nnr = size(nonrep,1);

ntrial = ntrial + nnr;

ndemas = 4 ;
nvalor = ceil((ntrial+ndemas)^(1/3));
ncomb  = (nvalor)^3 ;
colores =ones(ncomb,3)  ;
Val = 0:1/(nvalor-1):1;
% Vectors n2,n,1
V0 = repmat(Val,1,nvalor*nvalor) ;
Vk =[] ;
for i = 1:nvalor
    V1 = Val(i)*ones(size(Val));
    V1 = [Vk V1];
    Vk = V1 ;
end
V1 = repmat(Vk,1,nvalor);
Vk = [];
for i = 1:nvalor
    V2 = Val(i)*ones(1,(nvalor*nvalor));
    V2 = [Vk V2];
    Vk = V2 ;
end

colores(:,1) = V2';
colores(:,2) = V1';
colores(:,3) = V0';


colores = colores(1:ncomb-1,:);
switch ALEATORIO
    case 'YES'
        per_n = randperm(ncomb-1);
        colores = colores(per_n,:);
end

colores = colores((1:ntrial),:);
%colores = colores(randperm(ntrial),:);

%%%%5
newcolores = zeros(size(colores));
icolornew = 0 ;
if isempty(nonrep)
else
    for icolor = 1:size(colores,1)
        irep = 1 ;
        REP = 0 ;
        while irep < nnr

            % NO yellow
            if nonrep(irep,:) == colores(icolor,:)
                %'REPEATED !!!!'
                REP = 1;
                break
            else
                irep = irep + 1;
                % non-repeated

            end

        end
        if REP ==  0
            icolornew = icolornew + 1;
            newcolores(icolornew,:) = colores(icolor,:);

        end
    end
    if icolornew < ntrialor
        disp('WARNING in ColoresMatrix --> Repeated colors')
        colores = [newcolores(1:icolornew,:);yellow] ;
    else
        colores = newcolores(1:icolornew,:) ;
    end
end


%%%%%



%%% Randomly --->

%%%%
