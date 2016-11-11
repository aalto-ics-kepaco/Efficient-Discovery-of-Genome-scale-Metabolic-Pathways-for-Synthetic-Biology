function calculateVariationDistanceProduct()

%comparison based on all solutions identified by efms and k best. The
%number of solutions differ.


%% for the huge matrix
load('hugematrixResults.mat');
network=load('hugematrix.mat');
name=fieldnames(network);
inputstruct=network.(name{1});
[rkeep, outputstruct, D, BCD]=SimplifyStructure4Calc(inputstruct);
distanceProduct=outputstruct.rxndists;
distanceProduct=distanceProduct(rkeep);
   

%% for the yeast network 
load('fixedYeastResult.mat');
network=load(strcat(struct,'.mat'));
name=fieldnames(network);
outputstruct=network.(name{1});
rkeep=find((outputstruct.rxnClassification==0));
BCD=outputstruct.BCD;


% directions=sign(outputstruct.vL+outputstruct.vU);
% %directions=directions(rkeep);
% %dirs=directions;
% toRemove=[find(outputstruct.rxnClassification==1); find(outputstruct.rxnClassification==2)];
% directions(toRemove)=[];
% Rsize=vec2mat(outputstruct.distuni,size(outputstruct.R,2));
% Rsize(:,toRemove)=[];
% 
% unidirectionalSolutions=zeros(size(Solutions,1),sum(directions~=0)+2*sum(directions==0));
% distuni=zeros(1,sum(directions~=0)+2*sum(directions==0));
% k=1;
% BCDN=zeros(size(BCD,1),sum(directions~=0)+2*sum(directions==0));
% for i=1:size(BCD,2)
%     unidirectionalSolutions(:,k)=Solutions(:,i); 
%     BCDN(:,k)=BCD(:,i);
%     if (directions(i)==1)
%         distuni(k)=Rsize(1,i);
%          k=k+1;
%     elseif (directions(i)==-1)
%         distuni(k)=Rsize(2,i);
%         k=k+1;
%     else 
%         k=k+1;
%         distuni(k-1)=Rsize(1,i);
%         distuni(k)=Rsize(2,i);
%         unidirectionalSolutions(:,k)=-1*Solutions(:,i);
%         BCDN(:,k)=-1*BCD(:,i);
%         k=k+1;
%     end
% end


%% both cases 
% totalDist1=sum(outputstruct.rxndists(rkeep)==1);
% dist1ReactionsOld=find(outputstruct.rxndists(rkeep)==1);
% 

 dist1Reactions=find(distanceProduct==1);
 dist2Reactions=find(distanceProduct==2);
 dist3Reactions=find(distanceProduct==3);
 dist4Reactions=find(distanceProduct==4);
% dist1Reactions=find(distuni==1);
% dist2Reactions=find(distuni==2);
% dist3Reactions=find(distuni==3);
% dist4Reactions=find(distuni==4);
% totalDist2=sum(outputstruct.rxndists(rkeep)==2);
%dist2ReactionsOld=find(outputstruct.rxndists(rkeep)==2);
% 
% 
% totalDist3=sum(outputstruct.rxndists(rkeep)==3);
%% 
% totalDist4=sum(outputstruct.rxndists(rkeep)==4);
 %dist4ReactionsOld=find(outputstruct.rxndists(rkeep)==4);

%binaries=abs(unidirectionalSolutions)>1e-6;
binaries=abs(Solutions)>1e-6;
%binariesOld=abs(Solutions)>1e-6;
% solsDist1=nnz(sum(binariesOld(:,dist1Reactions)))/totalDist1;
% solsDist2=nnz(sum(binaries(:,dist2Reactions)))/totalDist2;
% solsDist3=nnz(sum(binaries(:,dist3Reactions)))/totalDist3;
% solsDist4=nnz(sum(binaries(:,dist4Reactions)))/totalDist4;
% solsDist=[solsDist1; solsDist2; solsDist3; solsDist4]; 


solsDist1=size(unique(binaries(:,dist1Reactions),'rows'),1);
%solsDist2=size(unique(binaries(:,dist2Reactions),'rows'),1)/totalDist2;
%solsDist3=size(unique(binaries(:,dist3Reactions),'rows'),1)/totalDist3;
%solsDist4=size(unique(binaries(:,dist4Reactions),'rows'),1)/totalDist4;

%% yeast 
%load('/home/milievsk/Documents/copenhagen /MATLAB/Comparison_EFMs/suppl_data/SuccinateCase/productConstraint/yeast/fixedYeastallVorigSize.mat');


%% huge
%load('/home/milievsk/Documents/copenhagen /MATLAB/Comparison_EFMs/suppl_data/SuccinateCase/productConstraint/hugematrixallVorigSize.mat');
load('/home/milievsk/Documents/copenhagen /MATLAB/Comparison_EFMs/suppl_data/SuccinateCase/productConstraint/yeast/fixedYeastallVorigSize.mat');

allVorigSizeBin=allVorigSize>1e-6;

%efmsDist1=nnz(sum(allVorigSizeBin(:,dist1ReactionsOld)))/totalDist1;
%efmsDist2=nnz(sum(allVorigSizeBin(:,dist2ReactionsOld)))/totalDist2;
%efmsDist3=nnz(sum(allVorigSizeBin(:,dist3ReactionsOld)))/totalDist3;
%efmsDist4=nnz(sum(allVorigSizeBin(:,dist4ReactionsOld)))/totalDist4;
%efmsDist=[efmsDist1; efmsDist2; efmsDist3; efmsDist4]; 
% the number of efms not sqrt size 
sqrtsize=size(allVorigSize,1);

% solsDist1sqrt=nnz(sum(binaries(1:sqrtsize,dist1Reactions)))/totalDist1;
% solsDist2sqrt=nnz(sum(binaries(1:sqrtsize,dist2Reactions)))/totalDist2;
% solsDist3sqrt=nnz(sum(binaries(1:sqrtsize,dist3Reactions)))/totalDist3;
% solsDist4sqrt=nnz(sum(binaries(1:sqrtsize,dist4Reactions)))/totalDist4;
% solsSqrt=[solsDist1sqrt; solsDist2sqrt; solsDist3sqrt; solsDist4sqrt]
%compute the portion of patways distance 1+2

reactions1plus2=[dist1Reactions' ;dist2Reactions'];
reactions1plus2Old=[dist1ReactionsOld ;dist2ReactionsOld];

% there are solutions with more than 1 dist 1 reaction
sols1Pos=sum(binaries(:,dist1Reactions),2)==1;
sols1Possqrt=sum(binaries(1:sqrtsize,dist1Reactions),2)==1;
%sols1=nnz(sum(binaries(:,dist1Reactions)));
%sols1sqrt=nnz(sum(binaries(1:sqrtsize,dist1Reactions)));
sols1=size(unique(binaries(sols1Pos,dist1Reactions),'rows'),1);
sols1sqrt=size(unique(binaries(sols1Possqrt,dist1Reactions),'rows'),1)

sols1plus2 = size(unique(binaries(:,reactions1plus2),'rows'),1)
sols1plus2sqrt = size(unique(binaries(1:sqrtsize,reactions1plus2),'rows'),1)
%efms1=size(unique(allVorigSizeBin(:,dist1ReactionsOld),'rows'),1)
efms1plus2= size(unique(allVorigSizeBin(:,reactions1plus2Old),'rows'),1)
efms1=nnz(sum(allVorigSizeBin(:,dist1ReactionsOld)));

reactions123=[reactions1plus2; dist3Reactions'];
reactions123Old=[reactions1plus2Old; dist3ReactionsOld];


sols123=size(unique(binaries(:,reactions123),'rows'),1)
sols123sqrt=size(unique(binaries(1:sqrtsize,reactions123),'rows'),1)

efms123= size(unique(allVorigSizeBin(:,reactions123Old),'rows'),1)

reactionsAll=[reactions123; dist4Reactions'];
reactionsAllOld=[reactions123Old; dist4ReactionsOld];

sols1234=size(unique(binaries(:,reactionsAll),'rows'),1)
sols1234sqrt=size(unique(binaries(1:sqrtsize,reactionsAll),'rows'),1)

efms1234=size(unique(allVorigSizeBin(:,reactionsAllOld),'rows'),1)

end