function  enumerateKdistance4SuccCaseYeast(struct)
%enumerate the first K flux solutions, distance 4 from product 


%%load('succinate.mat');
%load('majanetworks.mat',product);
%network=load('allnetworks.mat',product)
%product_name=fieldnames(network);
%struct=network.(product_name{1});

% add path to cplex matlab libraries 

addpath('/scratch/work/ilievsm1/ibm/cplex/matlab/x86-64_linux/');
addpath('/scratch/work/ilievsm1/ibm/cplex/examples/src/matlab/');
network=load(strcat(struct,'.mat'));
name=fieldnames(network);
outputstruct=network.(name{1});
%[rkeep, outputstruct, D, BCD]=SimplifyStructure4Calc(inputstruct);

%take BCD as the reaction matrix
react_mat=outputstruct.BCD;
%change directions between SR and CR
directions=sign(outputstruct.vL+outputstruct.vU);
directions=directions(outputstruct.rxnClassification==0);
%change to only forward and rev reactions
dirs=directions;
react_mat(:,directions<0)=-1*react_mat(:,directions<0);
directions(directions<0)=1;
rkeep=find((outputstruct.rxnClassification==0));

%read in constraints (not for yeast net)
%constr=outputstruct.shadowClassification;
%constr(outputstruct.shadowClassification==-2)=[];

%for the yeast net
constr=outputstruct.compoundClassification;

constraints=cell((size(constr)));
for i=1:size(constr,1)
    switch (constr(i))
        
        case 1
            constraints(i)=cellstr('substrate');
        case 0
            constraints(i)=cellstr('intermediate');
        case 2 
            constraints(i)=cellstr('product');
        %case -2
         %   constraints(i)=cellstr('cofactor');
    end
end

%add absolute vars and indicator variables
react_mat_indicator=[react_mat zeros(size(react_mat)) zeros(size(react_mat))];
[rows,columns]=size(react_mat);

orig_vars=linspace(1,columns,columns);
abs_vars=linspace(columns+1,2*columns,columns);
bins=linspace(2*columns+1,size(react_mat_indicator,2),columns);

% no cofactors
substrate_size=sum(strncmp((constraints),'substrate',length('substrate')));
product_size=sum(strncmp((constraints),'product',length('product')));
ineq_s=substrate_size+product_size;

% for yeast sparse netowrk use sprank (D) instead od rank(D).
nullity=size(outputstruct.D,2)-rank(outputstruct.D);

%oneineq constraint for substrate, 2*columns for abs and 2*columns for ind.
A=zeros(ineq_s+1+2*columns+2*columns,size(react_mat_indicator,2));
%equalities for intermediates 
Aeq=zeros(rows-ineq_s,size(react_mat_indicator,2));
k=1;
keq=1;
b=zeros(1,size(A,1));
beq=zeros(1,size(Aeq,1));

%sum of substrates
sum_s=zeros(1,size(react_mat,2));
for i=1:size(react_mat,1)

   % if (ismember(cellstr('cofactor'),constraints(i)))
    %    A(k,:)=[-1*react_mat(i,:) zeros(1,columns) zeros(1,columns)];
     %   b(k)=20;
      %  k=k+1;
       
    %end
    
    if (ismember(cellstr('substrate'),constraints(i)))
        A(k,:)=[react_mat(i,:) zeros(1,columns) zeros(1,columns)];
        b(k)=0;
        k=k+1;
         %sum of substrates
        sum_s=sum_s+react_mat(i,:);
    end
    
    if (ismember(cellstr('product'),constraints(i)))
        A(k,:)=[-1*react_mat(i,:) zeros(1,columns) zeros(1,columns)];
        b(k)=-1;
        k=k+1;
    end
    
    if ((ismember(cellstr('intermediate'),constraints(i))))
        Aeq(keq,:)=[react_mat(i,:) zeros(1,columns) zeros(1,columns)];
        beq(keq)=0;
        keq=keq+1;
    end
   
end

%sum of substrates constraint
A(k,:)=[sum_s zeros(1,size(react_mat,2)) zeros(1,columns)];
b(k)=-1; 
k=k+1;

%abs inequalities
for j=1:columns
    A(k,[j columns+j])=[-1 -1];
    A(k+1,[j columns+j])=[1 -1];
    b([k k+1])=[0 0];
    k=k+2;
end

%indicator constraints
for j=1:columns
    
    A(k,[columns+j 2*columns+j])=[1 -30];
    A(k+1,[columns+j 2*columns+j])=[-1 1];
    %%%
    b([k k+1])=[0 0];
    k=k+2;
    
end

lb(orig_vars)=-inf; lb(bins)=0; lb(abs_vars)=0;
lb(orig_vars(directions==1))=0;

ub(orig_vars)=inf; ub(abs_vars)=inf; ub(bins)=1;

%cost
cost=zeros(1,size(react_mat_indicator,2));
cost(columns+1:2*columns)=ones(1,columns);

%modifications for cplex
  f=cost';
  Aineq=A;
  bineq=b';
  beq=beq';
  ctype(orig_vars)='C'; ctype(abs_vars)='C'; ctype(bins)='B';
  
  running_times=zeros(nullity,1);
  tic;
  [x,~,~,~]=cplexmilp(f,Aineq,bineq,Aeq,beq,[],[],[],lb,ub,ctype);
  running_times(1)=toc;
  % create and mps model 
  %BuildMPS(A, b, Aeq, beq, cost,lb, ub, 'SumFluxes', 'Binary',bins,'MPSfilename',strcat(struct,'.1mps'));
   
   
  Solutions=zeros(nullity,columns);
  SolutionsBinary=zeros(nullity,columns);
  Objectives=zeros(1,nullity);
  Solutions(1,:)=x(1:columns);
  SolutionsBinary(1,:)=x(2*columns+1:end);
  Objectives(1)=sum(x(columns+1:2*columns));
  %running_times
  
  %calculate k=nullity solutions 
 for i=2:nullity
    
    A(k,:)=zeros(1,size(A,2));
    distanceProduct=outputstruct.rxndists;
    distanceProduct=distanceProduct(rkeep);
   
    % change between SR and CR   
    indicators=x(2*columns+1:end)';
    indicators(distanceProduct>4)=0;
    A(k,(2*columns+1:end))=indicators;
    
    %change for different exclusion constraint 
    b(k)=sum(indicators)-1;
    %b(k)=sum(indicators)-2;
    %b(k)=sum(indicators)-3;
%   
    k=k+1;
    Aineq=A;
    bineq=b';
    tic;
   [x,~,~,~]=cplexmilp(f,Aineq,bineq,Aeq,beq,[],[],[],lb,[],ctype);
   running_times(i)=toc
   BuildMPS(A, b, Aeq, beq, cost,lb, ub, 'SumFluxes', 'Binary',bins,'MPSfilename',strcat(struct,num2str(i),'.mps'));

    Solutions(i,:)=x(1:columns);
    SolutionsBinary(i,:)=x(2*columns+1:end);
    %fprintf('%d',i);
    Objectives(i)=sum(x(columns+1:2*columns));
end

%change fluxes for directions 
for i=1:size(directions,1)
    if (dirs(i)==-1)
        Solutions(:,i)=-1*Solutions(:,i);
    end
end
%time=toc;
save(strcat(struct,'Results'),'running_times','Solutions','Objectives','SolutionsBinary');
%save('succinatenoCofactorsAbs','time','Solutions','Objectives');
%save('glycolatenoCofactorsAbs','time','Solutions','Objectives');
%change react_mat directions 

%post-processing, presenting the results
%top50Solutions=Solutions(1:50,:);
%activePositionsTop20=find(sum(abs(top50Solutions)));

%change between SR and CR
%activeDist4=intersect(find(distanceProduct<=4),activePositionsTop20);
%activeDist4=intersect(find(distanceProduct<=4),activePositionsTop20);
%activeDistances=[distanceProduct(activeDist4)' ; int8(top20Solutions(:,activeDist4))];
%activeDistances=[distanceProductExchange(activeDist4)' ; (top50Solutions(:,activeDist4))];
%activeReactions=reactionIDs(activeDist4);

%[a,b]=sortrows(activeDistances',1);
%activeDistances=a'; activeReactions=activeReactions(b);
%save(strcat(inputstruct,'activeReactions'),'activeDistances','activeReactions');
%net=react_mat*top20Solutions';
