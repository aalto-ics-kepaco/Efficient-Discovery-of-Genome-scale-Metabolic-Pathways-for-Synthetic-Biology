function resp = classifyrxndir(vL, vU)
global myEpsilon;
resp = cell(size(vL));
resp((vL>-myEpsilon) & (vU<myEpsilon)) = {'closed'};
resp((vL<-myEpsilon) & (vU<myEpsilon)) = {'backward'};
resp((vL>-myEpsilon) & (vU>myEpsilon)) = {'forward'};
resp((vL<-myEpsilon) & (vU>myEpsilon)) = {'reversible'};
end