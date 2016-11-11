function displayinfoaboutstruct(inputstruct)

global myEpsilon;
myEpsilon=0.1;

%% populate local variables
% shadowKEGGCID = inputstruct.shadowKEGGCID; 
% shadowCIDnames = inputstruct.shadowCIDnames; 
% listofKEGGCID = inputstruct.listofKEGGCID; 
% listofCIDnames = inputstruct.listofCIDnames;
compoundClassification = inputstruct.compoundClassification; 
shadowClassification = inputstruct.shadowClassification;
R = inputstruct.R;
SR = inputstruct.SR;
CR = inputstruct.CR;
vL = inputstruct.vL;
vU = inputstruct.vU;
rxnClassification = inputstruct.rxnClassification; 

%% output info
[nspecies,nreactions] = size(CR);
disp(['CR has ',mat2str(nspecies),' metabolites and ',mat2str(nreactions),' reactions']);
[nspecies,nreactions] = size(SR);
disp(['SR has ',mat2str(nspecies),' metabolites and ',mat2str(nreactions),' reactions']);
[nspecies,nreactions] = size(R);
disp(['R has ',mat2str(nspecies),' metabolites and ',mat2str(nreactions),' reactions']);
substrates = compoundClassification==1;
products = compoundClassification==2;
cofactors = shadowClassification==-2;
intermediates = compoundClassification==0;
disp(['struct has ',mat2str(sum(substrates)),' substrates, ',mat2str(sum(products)),' products, ',mat2str(sum(cofactors)),' cofactors, and ',mat2str(sum(intermediates)),' intermediates']);
% disp('The substrates are')
% disp([listofKEGGCID(substrates),listofCIDnames(substrates)])
% disp('The cofactors are')
% disp([shadowKEGGCID(cofactors),shadowCIDnames(cofactors)])
% disp('The products are')
% disp([listofKEGGCID(products),listofCIDnames(products)]);
rxnclosed=((vL>-myEpsilon) & (vU<myEpsilon));
rxnbackward=((vL<-myEpsilon) & (vU<myEpsilon));
rxnforward=((vL>-myEpsilon) & (vU>myEpsilon));
rxnreversible=((vL<-myEpsilon) & (vU>myEpsilon));
disp(['struct has ',mat2str(sum(rxnreversible)),' reversible, ',mat2str(sum(rxnforward)),' forward, ',mat2str(sum(rxnbackward)),' backward, and ',mat2str(sum(rxnclosed)),' closed reactions']);

vexchangerxn = sum(abs(CR)>myEpsilon,1)==1;
vemptyrxn = sum(abs(CR)>myEpsilon,1)==0;
hasprodex = sum(abs(CR(shadowClassification==2,vexchangerxn))>myEpsilon,2)>0;
hassubex = sum(abs(CR(shadowClassification==1,vexchangerxn))>myEpsilon,2)>0;
hascofex = sum(abs(CR(shadowClassification==-2,vexchangerxn))>myEpsilon,2)>0;
disp(['CR has ',mat2str(sum(vemptyrxn)),' empty and ',mat2str(sum(vexchangerxn)),' exchange reactions']);
disp(['CR has ',mat2str(sum(hassubex)),' subex, ',mat2str(sum(hasprodex)),' prodex, and ',mat2str(sum(hascofex)),' cofex reactions']);
rpot = (vemptyrxn | vexchangerxn);
rpot = rpot(:,rxnClassification>-1);

vexchangerxn = sum(abs(SR)>myEpsilon,1)==1;
vemptyrxn = sum(abs(SR)>myEpsilon,1)==0;
hasprodex = sum(abs(SR(shadowClassification==2,vexchangerxn))>myEpsilon,2)>0;
hassubex = sum(abs(SR(shadowClassification==1,vexchangerxn))>myEpsilon,2)>0;
hascofex = sum(abs(SR(shadowClassification==-2,vexchangerxn))>myEpsilon,2)>0;
disp(['SR has ',mat2str(sum(vemptyrxn)),' empty and ',mat2str(sum(vexchangerxn)),' exchange reactions']);
disp(['SR has ',mat2str(sum(hassubex)),' subex, ',mat2str(sum(hasprodex)),' prodex, and ',mat2str(sum(hascofex)),' cofex reactions']);

vexchangerxn = sum(abs(R)>myEpsilon,1)==1 & rpot;
vemptyrxn = sum(abs(R)>myEpsilon,1)==0 & rpot;
hasprodex = sum(abs(R(compoundClassification==2,vexchangerxn))>myEpsilon,2)>0;
hassubex = sum(abs(R(compoundClassification==1,vexchangerxn))>myEpsilon,2)>0;
hascofex = sum(abs(R(compoundClassification==-2,vexchangerxn))>myEpsilon,2)>0;
disp(['R has ',mat2str(sum(vemptyrxn)),' empty and ',mat2str(sum(vexchangerxn)),' exchange reactions']);
disp(['R has ',mat2str(sum(hassubex)),' subex, ',mat2str(sum(hasprodex)),' prodex, and ',mat2str(sum(hascofex)),' cofex reactions']);


% prodsel = rxnClassification==2;
% substratesel = rxnClassification==1;
% disp(['There are ',mat2str(sum(substratesel)),' substrate reactions and ',mat2str(sum(prodsel)),' product reactions']);

end