function [rkeep, outputstruct,D,BCD] = SimplifyStructure4Calc(inputstruct)

%% change alternate products to intermediates
%     'C00042'    'Succinate'  
selprod = inputstruct.shadowClassification==2;
selsucc = ismember(inputstruct.shadowKEGGCID,{'C00042'});
shadowClassification = inputstruct.shadowClassification.*1;
shadowClassification(selprod & ~selsucc)=0;

%% make complete structure object
modif.shadowKEGGCID = inputstruct.shadowKEGGCID;
modif.shadowCIDnames = inputstruct.shadowCIDnames;
modif.shadowClassification = shadowClassification;
modif.SR = inputstruct.SR;
modif.vL = inputstruct.vL;
modif.vU = inputstruct.vU;
modif.listofKEGGRID = inputstruct.listofKEGGRID;
%oldFolder = cd('C:\DATA\2015\Tessa\cases\');
outputstruct = recreateandsplit(modif);
%cd(oldFolder);

%% remove all exchange reactions
rkeep = inputstruct.rxnClassification==0;
BCD = outputstruct.R(:,rkeep).*1;
D = BCD(outputstruct.compoundClassification==0,:).*1;
outputstruct.BCD = BCD;
outputstruct.D = D;

end