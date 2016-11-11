function modif = recreateandsplit(loadX)
% input compatible with loadX and saveX
% output compatible with R-SR-CR
% this file partially replaces fixnetqorks.m and convert.m

%% load from struct
%compoundClassification = loadX.compoundClassification;
%listofCIDnames         = loadX.listofCIDnames;
%listofKEGGCID          = loadX.listofKEGGCID;
listofKEGGRID          = loadX.listofKEGGRID;
%listofParEq            = loadX.listofParEq;
vL                     = loadX.vL;
vU                     = loadX.vU;
shadowCIDnames         = loadX.shadowCIDnames;
shadowClassification   = loadX.shadowClassification;
shadowKEGGCID          = loadX.shadowKEGGCID;
%R                      = loadX.R;
SR                     = loadX.SR.*1;

%% recreate and split
% generate appropriate numerical maximum
global myInfinity
myInfinity = 300;
% generate appropriate zero tolerance
global myEpsilon
myEpsilon = 1e-6;

% %% analyze SR
% vexchangerxn = sum(abs(SR)>myEpsilon,1)==1;
% vemptyrxn = sum(abs(SR)>myEpsilon,1)==0;
% nemptyrxn = sum(vemptyrxn);
% nexchangerxn = sum(vexchangerxn);
% ncof = sum(shadowClassification==-2);
% nsub = sum(shadowClassification==1);
% nprod = sum(shadowClassification==2);
% hasprodex = sum(abs(SR(shadowClassification==2,vexchangerxn))>myEpsilon)>0;
% hassubex = sum(abs(SR(shadowClassification==1,vexchangerxn))>myEpsilon,2)>0;
% hascofex = sum(abs(SR(shadowClassification==-2,vexchangerxn))>myEpsilon,2)>0;

% %% remove closed and empty reactions
% rclosed = (vL>-myEpsilon & vU<myEpsilon);
% SR(:,rclosed)=0; % this causes it to be removed below

%% remove exchange reactions
rxnkeep = sum(abs(SR)>myEpsilon,1)>1;
CR = SR(:,rxnkeep);
vL = vL(rxnkeep');
vU = vU(rxnkeep');
listofKEGGRID = listofKEGGRID(rxnkeep');
clear rxnskeep SR

% %% remove species not used in any reaction
% skeep = sum(abs(CR)>myEpsilon,2)>0;
% CR = CR(skeep,:);
% shadowCIDnames         = shadowCIDnames(skeep);
% shadowClassification   = shadowClassification(skeep);
% shadowKEGGCID          = shadowKEGGCID(skeep);
nspcs = size(CR,1);

%% add exchange reactions
cofex = [];
lstcofex = {};
rxnClassification = zeros(size(vL));
for k = 1:nspcs
    if(shadowClassification(k)~=0)
        rxne = zeros(nspcs,1);
        rxne(k) = -1;
        if(shadowClassification(k)==-2)
            cofex = [cofex.*1,rxne.*1];
            lstcofex = [lstcofex;strcat('x',shadowKEGGCID(k))];
        else
            CR = [CR.*1,rxne.*1];
            if(shadowClassification(k)==2)
                vL = [vL;0];
                vU = [vU;myInfinity];
                rxnClassification=[rxnClassification;2];
            else %'==1 assumed
                vL = [vL;-myInfinity];
                vU = [vU;0];%myInfinity
                rxnClassification=[rxnClassification;1];
            end
            listofKEGGRID = [listofKEGGRID;strcat('x',shadowKEGGCID(k))];
        end
    end
end
CR = [CR.*1,cofex.*1];
allKEGGRID = [listofKEGGRID;lstcofex];
sones=ones(size(cofex,2),1);
avL = [vL;sones.*-myInfinity];
avU = [vU;sones.*myInfinity];
allrxnClassification=[rxnClassification.*1;sones.*-2];
clear nspcs rxne k cofex lstcofex sones

%% split network
metkeep = (shadowClassification>-myEpsilon);
rxnkeep = (allrxnClassification>-myEpsilon);
% disp([sum(metkeep),sum(~metkeep),sum(rxnkeep),sum(~rxnkeep)])
R = CR(metkeep,rxnkeep).*1;
listofKEGGCID = shadowKEGGCID(metkeep);
listofCIDnames = shadowCIDnames(metkeep);
compoundClassification = shadowClassification(metkeep);

tempR = CR.*1;
tempR(:,allrxnClassification~=0)=0;
SR = tempR(:,rxnkeep).*1;

metkeep = (shadowClassification~=0);
S = tempR(metkeep,:).*1;
netKEGGCID = shadowKEGGCID(metkeep);
netCIDnames = shadowCIDnames(metkeep);
netClassification = shadowClassification(metkeep);
clear metkeep tempR

[rxndists, distuni] = calcdists(R,vL.*1,vU.*1,compoundClassification.*1);

%% reaction directions
rxndir = zeros(size(vL));
rxndir((vL<-myEpsilon) & (vU<myEpsilon)) = -1;
rxndir((vL>-myEpsilon) & (vU>myEpsilon)) = 1;
% disp(['The full network has ',mat2str(size(CR,1)),' metabolites and ',mat2str(size(CR,2)),' reactions']);
% disp(['There are ',mat2str(sum(shadowClassification==1)),' substrates, ',mat2str(sum(shadowClassification==2)),' products, and ',mat2str(sum(shadowClassification==-2)),' cofactors']);
% disp(['There are ',mat2str(sum(rxndir==0)),' reversible and ',mat2str(sum(rxndir~=0)),' irreversible reactions']);
% %disp(['The product is ',shadowKEGGCID(shadowClassification==2),shadowCIDnames(shadowClassification==2)]);

%% save to struct
modif.allKEGGRID = allKEGGRID;
modif.CR = CR;
modif.avL = avL;
modif.avU = avU;
modif.allrxnClassification = allrxnClassification;

modif.shadowKEGGCID = shadowKEGGCID;
modif.shadowCIDnames = shadowCIDnames;
modif.shadowClassification = shadowClassification;
modif.SR = SR;
modif.vL = vL;
modif.vU = vU;
modif.rxndir = rxndir;
modif.rxnClassification = rxnClassification;
modif.listofKEGGRID = listofKEGGRID;
modif.R = R;
modif.listofKEGGCID = listofKEGGCID;
modif.listofCIDnames = listofCIDnames;
modif.compoundClassification = compoundClassification;
modif.S = S;
modif.netKEGGCID = netKEGGCID;
modif.netCIDnames = netCIDnames;
modif.netClassification = netClassification;
modif.rxndists = rxndists;
modif.distuni = distuni;

end
