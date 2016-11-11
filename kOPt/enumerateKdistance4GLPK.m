function  enumerateKdistance4GLPK(struct)
%enumerate the first K, reaction space only, distance 4 from product 

% input:
% struct - char name of the reaction matrix, ex. 'candidate159x282'

network=load(strcat(struct,'.mat'));
name=fieldnames(network);
inputstruct=network.(name{1});
[rkeep, outputstruct, D, BCD]=SimplifyStructure4Calc(inputstruct);

%take BCD as the reaction matrix
react_mat=BCD;

%directions=sign(struct.vL+struct.vU);
directions=sign(outputstruct.vL+outputstruct.vU);
directions=directions(rkeep);
%change to only forward and rev reactions, modify the reaact_matrix as well
dirs=directions;
react_mat(:,directions<0)=-1*react_mat(:,directions<0);
directions(directions<0)=1;

%% read in met. classification and create constraints

constr=outputstruct.shadowClassification;
constr(outputstruct.shadowClassification==-2)=[];
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% add absolute vars and indicator variables to the reaction matrix 
% create constraints for the type of variables 


react_mat_indicator=[react_mat zeros(size(react_mat)) zeros(size(react_mat))];
[rows,columns]=size(react_mat);

% orig_vars refers to the reaction fluxes
orig_vars=linspace(1,columns,columns);
abs_vars=linspace(columns+1,2*columns,columns);
bins=linspace(2*columns+1,size(react_mat_indicator,2),columns);

%% build the MILP model 
% no cofactors 
substrate_size=sum(strncmp((constraints),'substrate',length('substrate')));
product_size=sum(strncmp((constraints),'product',length('product')));

ineq_s=substrate_size+product_size;
nullity=size(D,2)-rank(D);

%oneineq constraint for substrate, 2*columns for abs and 2*columns for
%indicator vars 
A=zeros(ineq_s+1+2*columns+2*columns,size(react_mat_indicator,2));
%equalities for intermediates 
Aeq=zeros(rows-ineq_s,size(react_mat_indicator,2));
k=1;
keq=1;
b=zeros(1,size(A,1));
beq=zeros(1,size(Aeq,1));

%sum of substrates constraint 
sum_s=zeros(1,size(react_mat,2));
for i=1:size(react_mat,1)
    
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
 % lower and upper boun and cost 
lb(orig_vars)=-inf; lb(bins)=0; lb(abs_vars)=0;
lb(orig_vars(directions==1))=0;
ub(orig_vars)=inf; ub(abs_vars)=inf; ub(bins)=1;
cost=zeros(1,size(react_mat_indicator,2));
cost(columns+1:2*columns)=ones(1,columns);

%% modifications for cplex, running the optimization 
  f=cost';
  Aineq=A;
  bineq=b';
  beq=beq';
  ctype(orig_vars)='C'; ctype(abs_vars)='C'; ctype(bins)='B';
  
  %% modification for GLPK 
  cGLPK=f;
  aGLPK=[Aeq; Aineq];
  kGLPK=k+keq-1;
  bGLPK=[beq ; bineq];
  vartypeGLPK=ctype;
  ctypeGLPK=[repmat('S',size(Aeq,1),1);  repmat('U',size(Aineq,1),1) ];
  sense=1;
  param.msglev=3;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  running_times=zeros(nullity,1);
  tic;
  %[x,~,~,~]=cplexmilp(f,Aineq,bineq,Aeq,beq,[],[],[],lb,ub,ctype);
  [xmin,fmin,~, ~]=glpkcc(cGLPK,aGLPK,bGLPK,lb,ub,ctypeGLPK',vartypeGLPK,sense,param);
  
  running_times(1)=toc;
  % create and mps model, optional; BuildMPS function required 
  %BuildMPS(A, b, Aeq, beq, cost,lb, ub, 'SumFluxes', 'Binary',bins,'MPSfilename',strcat(struct,'.1mps'));
    
  Solutions=zeros(nullity,columns);
  SolutionsBinary=zeros(nullity,columns);
  Objectives=zeros(1,nullity);
  Solutions(1,:)=xmin(1:columns);
  SolutionsBinary(1,:)=xmin(2*columns+1:end);
 % Objectives(1)=sum(xmin(columns+1:2*columns));
  Objectives(1)=fmin;
 %running_times
  
  %calculate k=nullity solutions 
  distanceProduct=outputstruct.rxndists;
  distanceProduct=distanceProduct(rkeep);
 for i=2:nullity
    
    aGLPK(kGLPK,:)=zeros(1,size(aGLPK,2));
   
    % adding the exclusion constraint 
    indicators=xmin(2*columns+1:end)';
    indicators(distanceProduct>4)=0;
    aGLPK(kGLPK,(2*columns+1:end))=indicators;
    
    %change for different exclusion constraint,  
    bGLPK(kGLPK)=sum(indicators)-1;
    ctypeGLPK(kGLPK)='U';
    %b(k)=sum(indicators)-2;
    %b(k)=sum(indicators)-3;
%   
    kGLPK=kGLPK+1;
   % Aineq=A;
    %bineq=b';
    tic;
    %[x,~,~,~]=cplexmilp(f,Aineq,bineq,Aeq,beq,[],[],[],lb,[],ctype);
     [xmin,fmin,~, ~]=glpkcc(cGLPK,aGLPK,bGLPK,lb,ub,ctypeGLPK',vartypeGLPK,sense,param);
    running_times(i)=toc;
   %BuildMPS(A, b, Aeq, beq, cost,lb, ub, 'SumFluxes', 'Binary',bins,'MPSfilename',strcat(struct,num2str(i),'.mps'));

    Solutions(i,:)=xmin(1:columns);
    SolutionsBinary(i,:)=xmin(2*columns+1:end);
    Objectives(i)=fmin;
end

%change fluxes for directions 
for i=1:size(directions,1)
    if (dirs(i)==-1)
        Solutions(:,i)=-1*Solutions(:,i);
    end
end
%time=toc;
save(strcat(struct,'ResultsGLPK'),'running_times','Solutions','Objectives','SolutionsBinary');
