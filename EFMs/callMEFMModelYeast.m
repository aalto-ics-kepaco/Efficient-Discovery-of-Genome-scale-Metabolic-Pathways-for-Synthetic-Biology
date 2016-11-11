function  callMEFMModelYeast(struct)
    global allKEFMS
    global allVorigSize
%input:
% product - name of the struct that contain the product S matrix
% K - number of EFMS that we want to calculate
% place the file majanetworks.mat in the same working directory, or include
% the path

addpath('/scratch/work/ilievsm1/ibm/cplex/matlab/x86-64_linux/');
addpath('/scratch/work/ilievsm1/ibm/cplex/examples/src/matlab/');

%load('succinate.mat');
%load('majanetworks.mat')
%network=load('allnetworks.mat',product)
%product_name=fieldnames(network);
%struct=network.(product_name{1});
network=load(strcat(struct,'.mat'));
name=fieldnames(network);
outputstruct=network.(name{1});

nullity=size(outputstruct.D,2)-rank(outputstruct.D);
react_mat=outputstruct.BCD;

%react_mat=struct.CR;
directions=sign(outputstruct.vL+outputstruct.vU);
directions=directions(outputstruct.rxnClassification==0);
dirs=directions;
react_mat(:,directions<0)=-1*react_mat(:,directions<0);
directions(directions<0)=1;
%rkeep=find((outputstruct.rxnClassification==0));

%read in constraints
%constr=outputstruct.shadowClassification;
%constr(outputstruct.shadowClassification==-2)=[];

constr=outputstruct.compoundClassification;
constraints=cell((size(constr)));
for i=1:size(constr,1)
    switch (constr(i))
        
        case  1
            constraints(i)=cellstr('substrate');
        case  0
            constraints(i)=cellstr('intermediate');
        case  2 
            constraints(i)=cellstr('product');
        case -2
            constraints(i)=cellstr('cofactor');
    end
end

%intermediates, now only substrates and product
nonNegativeOnly=zeros(size(react_mat,1),sum(directions==1)+2*sum(directions==0));
k=1;
for i=1:size(react_mat,2)
    nonNegativeOnly(:,k)=react_mat(:,i);
    k=k+1;
 if (directions (i)==0)    
     nonNegativeOnly(:,k)=-1*react_mat(:,i);
    k=k+1;
 end
end


% S is the stoichiometrix matrix required as input for function
% mEFMModel.m
% creating the variables for cplex model here..
% check cplex specifications for details 
outputstruct.shadowClassification=outputstruct.shadowClassification(outputstruct.shadowClassification~=-2);
S=nonNegativeOnly(outputstruct.shadowClassification==0,:);
subs=nonNegativeOnly(outputstruct.shadowClassification==1,:);
prod=-1*nonNegativeOnly(outputstruct.shadowClassification==2,:);
sumSubs=sum(subs);
constrIni.A=[subs; prod; sumSubs];
constrIni.lhs=repmat(-inf,1,size(constrIni.A,1));
constrIni.rhs=[zeros(1,size(subs,1)) -1 -1]; 
constrIni.lb=repmat(-inf,1,size(constrIni.A,1));
constrIni.up=inf(1,size(constrIni.A,1));
constrIni.lb(directions==1)=0;


%constraint r6 for enumeration
% the function mEFMModel requires the user to enter the enumeration
% constraint manualy

%we save all the EFMs to a struct, and we load them in the struct exists 
if exist('allKEFMS.mat','file')
    %load the binary EFMs 
    load('allKEFMS.mat');
    load('allVorigSize.mat');
else
    %otherwise we generate the struct for the first EFM
    allKEFMS=zeros(1,size(S,2));
    allVorigSize=zeros(1,size(react_mat,2));
end

% loop for calculating K EFMs,
%running_time=zeros(NumEFMS,1);  %BuildMPS(A, b, Aeq, beq, cost,lb, ub, 'SumFluxes', 'Binary',bins,'MPSfilename','SumFluxes4.mps');

% max time is 12 hours 
max_running_time=43200;
running_times=zeros(nullity,1);
total_time=0;
counter=1;
while total_time <=max_running_time
     if counter>nullity
         break;
     end
     [opt_time,EFM1,EFM1bin, ~, ~, ~]=mEFMModel(S,constrIni,counter);
running_times(counter)=opt_time;
total_time=total_time+opt_time;
%cplex.writeModel(strcat(int2str(K),'EFM.lp'));
%cplex.Model.rhs;
allKEFMS(counter,:)=EFM1bin';


%save the flux vector, with bi-directional reactions
 vOriginalSize=zeros(size(react_mat,2),1);
       k=1;
       for i=1:size(react_mat,2)
            vOriginalSize(i)=EFM1(k);
            k=k+1;
            if(directions(i)==0)
                vOriginalSize(i)=vOriginalSize(i)-EFM1(k);
                k=k+1;
            end
       end
       
allVorigSize(counter,:)=vOriginalSize';
counter=counter+1;
end

%modify the directions for reverse only reactions
for i=1:size(directions,1)
   if (dirs(i)==-1)
      allVorigSize(:,i)=-1*allVorigSize(:,i);
   end
end
save(strcat(struct,'allKEFMS.mat'),'allKEFMS');
save(strcat(struct,'allVorigSize.mat'),'allVorigSize');
%modify to present only the reactions within distance 4 from product 
%activePositionsTop20=find(sum(abs(allVorigSize)));

%no cofactors distance from product 
% distanceProduct=struct.rxndists;
% distanceProductExchange=100*ones(size(react_mat,2),1);
% nonExchangeReactions= struct.rxnClassification~=-2;
% distanceProductExchange(nonExchangeReactions)=distanceProduct;
% 
% activeDist4=intersect(find(distanceProductExchange<=4),activePositionsTop20);
% activeDistances=[distanceProductExchange(activeDist4)'; (allVorigSize(:,activeDist4))];
% activeReactions=struct.listofKEGGRID(activeDist4);
% [a,b]=sortrows(activeDistances',1);
% activeDistances=a'; activeReactions=activeReactions(b);
save(strcat(struct,'running_time.mat'),'running_times');
end