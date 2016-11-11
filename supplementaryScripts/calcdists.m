function [dist, distuni] = calcdists(R,vL,vU,compoundClassification)
myEpsilon=1e-6;

%% obtain R ignoring cofactors

%% Make unidirectional reaction matrix by duplication of reversible reactions
nreactions = size(vL,1);
selvL = (vL<-myEpsilon); % reversible and backward-only
selvU = (vU>myEpsilon); % reversible and forward-only
uniR=[R(:,selvU),-R(:,selvL)];
map=eye(2*nreactions);
map=map(:,[selvU;selvL]);
mapuniRtoR = map;

%% create list of remaining metabolites
spcvisited = zeros(size(uniR,1),1);

%% set product as "reached"
spcvisited(compoundClassification==2)=1;

%% iterate
[rxnvisited] = iterate(uniR, spcvisited);

%% postprocess
nums = mapuniRtoR * rxnvisited;
distuni = nums.*1;
dist = reshape(nums,nreactions,2);
% size(nums);
dist(dist==0)=inf;
dist = min(dist,[],2);
dist(dist==inf)=0;
  
end

function [rxnvisited] = iterate(R, spcvisited)
%% initialize
rxnvisited = zeros(size(R,2),1);
Rm = R.*1;
Rp = R.*-1;
Rm(Rm<0)=0; % one side of all reactions
Rp(Rp<0)=0; % other side of all reactions

for iter = 1:200
    disp(['iteration ',mat2str(iter)]);
    
    %% get list of all remaining reactions having at least one substrate in the list of reached species
    rxnnew = ((Rm' * spcvisited)>0);
    
    %% remove reactions already visited
    rxnnew(rxnvisited>0)=false;

    %% if no new reactions, exit for
    n = sum(rxnnew);
    if(n==0)
        %disp('no new reactions');
        break;
    else
        %disp([mat2str(n),' new reactions']);
    end
    
    %% flag reactions as reached
    rxnvisited(rxnnew)=iter;

    %% get list of products in those reactions
    spcprod = (sum(Rp(:,rxnnew'),2)>0);

    %% remove species already visited
    spcprod(spcvisited>0)=false;
    
    %% update list of remaining species by setting values to iter
    spcvisited(spcprod)=(iter+1);
    
    %% if no new species, exit for
    n = sum(spcprod);
    if(n==0)
        %disp('no new species');
        break;
    else
        %disp([mat2str(n),' new species']);
    end
end

end