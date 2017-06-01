
[BuildMPS](https://se.mathworks.com/matlabcentral/fileexchange/19618-mps-format-exporting-tool?focused=5153282&tab=function) is a function for converting linear programming problem (MATLAB matrices) to MPS format, developed by  Bruno Luong. 
The function is only used for representing the linear programs into more comprehensible MPS format. 

Two versions of the kOpt included: 

enumerateKdistance4GLPK.m - to be used with GLPK optimization software.

1) enumerateKdistance4SuccinateKEGG.m
2) enumerateKdistance4SuccinateYeast.m, both to be used with CPLEX. 
Due to some minor differences in representation of the structures, 1) and 2) are used depending on wherther the yeast network or KEGG subnetworks are required.




