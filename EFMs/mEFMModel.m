%function [cplex,v,z] =mEFMModel(S,constrIni,K)
function [optimization_time,EFM1,EFM1bin,cplex, var, constr] = mEFMModel(S,constrIni,K)
     global allKEFMS
   
% mEFMModel   generates the CPLEX model to enumerate K Elementary Flux Modes 
%             calculated in the ligth of the stoichiometric matrix S
%             satisfying constraints in constrIni.
% 
%    [cplex var constr] = mEFMModel(S,constrIni,K)
% 
%    cplex is the structure with the generated optimization model.
%    var are the positions of each variable in the model: var.x,v,z,d,g.
%    constr is an structure with the position corresponding to each family 
%    of constraints: 
%       constr.r1 --> Equation 15
%       constr.r2 --> Equation 16 
%       constr.r3 --> Equation 10b
%       constr.r4 --> Equation 1, 3 and 14
%       constr.r5 --> Equation 10a
%       constr.r6 --> Equation 13
% 
%    The CPLEX parameters can be subsequently tuned.
%    Recursive constraints preventing the calculation of previously retrieved 
%    solutions are manually included by the user in the positions within the 
%    array constr.r6.


tic;
M = 1000;

cplex=Cplex('p1');
cplex.Model.sense='minimize';
cplex.Param.output.clonelog.Cur = -1;


nConstra = size(constrIni.A,1);
%indicator variables
z = 1:size(S,2);
%total number of variables
v = z(end)+1:(z(end)+size(S,2));
nEnd = v(end);
%epsilons  
d = reshape(nEnd+(1:size(S,2)*(nConstra-1)),size(S,2),nConstra-1);
nEnd = d(end);
%deltas 
g = reshape(nEnd+(1:size(S,2)*(nConstra-1)),size(S,2),nConstra-1);
nEnd = g(end);
x = reshape(nEnd+(1:((size(S,1)+1)*(nConstra-1))),size(S,1)+1,nConstra-1);
nVar = x(end);


r2 = +1:length(d);
r3 = r2(end)+1:r2(end)+length(z);
r4 = r3(end)+1:r3(end)+size(S,1)+nConstra;
nEnd = r4(end);
r1 = reshape(nEnd+(1:size(S,2)*(nConstra-1)),size(S,2),nConstra-1);
nEnd = r1(end);
r5 = nEnd+1:nEnd+length(z);
%if K>0
if K>1
    %r6= r5(end)+1:r5(end)+K
    r6 = r5(end)+1:r5(end)+K-1;
    nCons = r6(end);
else
    nCons = r5(end);
end


nel = sum(sum(S~=0))*2*3*nConstra + length(z)*2*nConstra + length(v) + 1;
A = spalloc(nCons,nVar,nel);
nel = 1 + length(r2) + length(r3);
lhs = spalloc(nCons,1,nel);
rhs = spalloc(nCons,1,nel);
lb = spalloc(nVar,1,length(x)+3); lb(x) = -inf; 
ub = inf(nVar,1); ub(z)=1;
obj = spalloc(nVar,1,0);
%variables are integer
ctype(x) = 'C';ctype(v) = 'C';ctype(d) = 'C';
ctype(g) = 'C';ctype(z) = 'B';


Sp = [S' zeros(size(S,2),1)]; Sp(:,end) = constrIni.A(end,:)';
for i = 1: nConstra -1   
    A(r1(:,i),x(:,i)) = Sp;
    ind = sub2ind(size(A),r1(:,i)',g(:,i)');    A(ind) = 1;
    ind = sub2ind(size(A),r1(:,i)',d(:,i)');    A(ind) = -1;    
    lhs(r1(:,i)) = constrIni.A(i,:)';
    rhs(r1(:,i)) = constrIni.A(i,:)';
   
    ind = sub2ind(size(A),r2,g(:,i)');    A(ind) = 1;
    ind = sub2ind(size(A),r2,d(:,i)');    A(ind) = 1;    
end
ind = sub2ind(size(A),r2,z);    A(ind) = M;
lhs(r2) = -inf; rhs(r2) = M;

ind = sub2ind(size(A),r3,z);    A(ind) = -M;
ind = sub2ind(size(A),r3,v);    A(ind) = 1;
lhs(r3) = -inf; rhs(r3) = 0;

A(r4(1:size(S,1)),v) = S;
A(r4(size(S,1)+1:end),v) = constrIni.A;
lhs(r4(size(S,1)+1:end)) = constrIni.lhs;
rhs(r4(size(S,1)+1:end)) = constrIni.rhs;

ind = sub2ind(size(A),r5,z);    A(ind) = -1;
ind = sub2ind(size(A),r5,v);    A(ind) = 1;
lhs(r5) = 0; rhs(r5) = inf;
%add constraint r6 for enumeration 
if K>1
    %exclusion on the binary vars
    A(r6,z)=allKEFMS;
    rhs(r6)=sum(allKEFMS,2)-1 ;  
    lhs(r6)=-inf;
end


cplex.Model.A=A;    cplex.Model.obj=obj;
cplex.Model.lhs=lhs;    cplex.Model.rhs=rhs;
cplex.Model.ub=ub;  cplex.Model.lb=lb;
cplex.Model.ctype = ctype;

cplex.Model.obj(z) = 1;

var.x = x;var.v = v;var.z = z;var.d = d;var.g = g;
constr.r1 = r1;constr.r2 = r2;constr.r3 = r3;constr.r4 = r4;constr.r5 = r5;
%if K>0
if K>1
    constr.r6 = r6;
end

%%%%% first version of the cplex model is generated 
model_generation_time=toc;
tic;
solution=cplex.solve;
optimization_time=toc;
x=solution.x
EFM1=x(v); 
EFM1bin=x(z);
end

