function calculatevariationfromsolutions()
% Autrhor: Peter Blomberg, VTT
% 19.8.2016
% 22.8.2016 modified

% todo: chk if yeast76succinate is the same as modifiedyeast76succinate YES, they are
% todo: calc variation for Maja's solutions

% global myInfinity
% myInfinity = 300;
global myEpsilon
myEpsilon=1e-6;

% YEAST
load fixedYeastResult.mat
load fixedYeast.mat 
inputstruct = outputstruct; % from fixedYeast.mat
%vector2filter = Solutions; % from fixedYeastResults.mat
vector2filter = SolutionsBinary; % from fixedYeastResults.mat
% The data is in unidirectional format

% optionally limit to the first 23 flux solutions (rows)
vector2filter = vector2filter(1:23,:);

% %% KEGG
% load hugematrixResults.mat
% load hugematrix.mat
% inputstruct = hugematrix; % from hugematrix.mat
% %vector2filter = Solutions; % from hugematrixResults.mat
% vector2filter = SolutionsBinary; % from hugematrixResults.mat
%
% % optionally limit to the first 35 flux solutions (rows)
% vector2filter = vector2filter(1:35,:);

%% convert from binary to Peter's unidirectional

% % The data is in bidirectional format
% selpos = vector2filter.*1;
% selpos(vector2filter<myEpsilon)=0;
% selneg = vector2filter.*1;
% selneg(vector2filter>-myEpsilon)=0;
% vector2filter = [selpos,-selneg];


%% continue calculations

% quality check
% max(abs(Solutions(SolutionsBinary<0.5)))
% min(abs(Solutions(SolutionsBinary>0.5)))

% fix near binary to binary
vector2filter(abs(vector2filter)<myEpsilon)=0;
vector2filter(abs(vector2filter-1)<myEpsilon)=1;

% calculate the stats
[uniquefiltered, filteredonly, ufcount, fosize, forank, fox] = filterfunction(inputstruct, vector2filter, 1);
disp([ufcount, fosize, fox, forank])
[uniquefiltered, filteredonly, ufcount, fosize, forank, fox] = filterfunction(inputstruct, vector2filter, 2);
disp([ufcount, fosize, fox, forank])
[uniquefiltered, filteredonly, ufcount, fosize, forank, fox] = filterfunction(inputstruct, vector2filter, 3);
disp([ufcount, fosize, fox, forank])
[uniquefiltered, filteredonly, ufcount, fosize, forank, fox] = filterfunction(inputstruct, vector2filter, 4);
disp([ufcount, fosize, fox, forank])

end


function [uniquefiltered, filteredonly, ufcount, fosize, forank, fox] = filterfunction(inputstruct, vector2filter, filterdist)
rxndists = inputstruct.distuni;
%size(rxndists) %5298x1  4682x1
c4uni = zeros(size(rxndists));
c4uni(rxndists<(filterdist+0.4) & rxndists>0.4)=1;
%size(c4uni) %5298x1 4682x1
%sum(c4uni) %6   25  89  284    343 1331    2190    2673
rkeep = inputstruct.rxnClassification==0;
%size(rkeep) %2649x1 2341x1
%sum(rkeep) %2642    2271
filter4uni4BCD = [rkeep;rkeep];
filter4uni4BCD = filter4uni4BCD(c4uni>0.5);
%size(filter4uni4BCD) %6x1   25x1    89x1    284x1   343x1
%sum(filter4uni4BCD) %6  25  89  284 343    1327    2182    2664

% vector is 1208x2271
% filter is 2673

filteredonly = vector2filter(:,filter4uni4BCD);
%size(filteredonly) %1271x6  1271x25 1271x89 1271x284    1208x343

[fosize,fox]=size(filteredonly);
forank=rank(filteredonly);
if(forank==0)
    filteredonly
end

uniquefiltered = unique(filteredonly,'rows');
%size(uniquefiltered) %127x6 231x25  322x89  842x284 325x343    1204x4327
%1208x2182
ufcount=size(uniquefiltered,1);

end