%% function to generate the plots
function generatePlotRunningTimes()
%final version%%%%%%%%%%%

% add the path to the location storing the kOpt results structures
%add the path to the location storing the EFMs results structures
%addpath('/home/milievsk/Documents/copenhagen /MATLAB/pathConstructor/BuildMPS/BuildMPS/distance4Results/Succinate9case /');
% suffix for the EFMs results structs is 'running_time.mat'
%suffix for tje kOPt results is 'Results.mat'
%add the path to the additional functions located in the
%SupplementaryScripts folder 

filesEFMS=dir('*running_time.mat');
[sqrtN,sorted_reactionsEFMs,running_timesEFMs,idxs]=generateSqrtRunningTimes(filesEFMS);
%figure1 = figure;

% Create axes
%axes1 = axes('Parent',figure1);
%hold(axes1,'on');

% Create plot
%plot
%%%%%
fig=figure;
semilogy(sqrtN(idxs),running_timesEFMs(idxs)/60,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0 0 1],...
'MarkerSize',10,...
    'Marker','square',...
    'LineWidth',2,...
    'LineStyle','--',...
    'Color',[0 1 0]);


% Create xlabel
%xlabel('number of reactions','FontSize',11);
xlabel('Number of pathways')
% Create ylabel
ylabel('Computation time (minutes)','FontSize',11);
%box(axes1,'on');

legend('EFMs');
hold on;

%% add the running time for k-best 
%load('TimeVsNumreactKbest');
%filesKbest=dir('/home/milievsk/Documents/copenhagen /MATLAB/pathConstructor/BuildMPS/BuildMPS/distance4Results/Succinate9case /*Results.mat');
filesKbest=dir('*Results.mat');

[sqrtN,sorted_reactionsKbest,running_timesKbest, idxsKbest]=generateSqrtRunningTimes(filesKbest);
%semilogy(sorted_reactionsKbest,running_timesKbest(idxsKbest)/60,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0 0 1],...
semilogy(sqrtN(idxsKbest),running_timesKbest(idxsKbest)/60,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0 0 1],...
'MarkerSize',10,...
    'Marker','square',...
    'LineWidth',2,...
    'LineStyle','--',...
    'Color',[1 0 0]);
legend('EFMs','kOpt')
%legend('kOpt')
%legend('k-best')
%print tif figure required for publication
print(fig,'ComTimes','-dtiffn');
end

%% function generate times for sqrt solutions
function [sqrtN,sorted_reactions, times, idxs] = generateSqrtRunningTimes(files)
times=zeros(length(files),1);
num_reactions=zeros(length(files),1);
sqrtN=zeros(length(files),1);
for i=1:length(files)
    load(files(i).name);
    %running_times(running_times==0)=[];

    %retrieve the nullity of the network 
    struct=files(i).name;  
    switch files(1).name((end-4):end) 
        case 'e.mat'
        struct(end-(length('running_time.mat')-1):end)=[];
        otherwise
        struct(end-(length('Results.mat')-1):end)=[];
    end
    network=load(strcat(struct,'.mat'));
    name=fieldnames(network);
    inputstruct=network.(name{1});
    [~,~,D,~]=SimplifyStructure4Calc(inputstruct);
   % size(BCD)
    nullity=size(D,2)-rank(D);
    % nullN(i)=nullity;
    %keep only the sqrt running times 
    sqRootNullity=round(sqrt(nullity));
    sqrtN(i)=sqRootNullity;
    running_times(sqRootNullity+1:end)=[];

    %save the time 
    times(i)=sum(running_times);
    num_reactions(i)=size(D,2);
end
[sorted_reactions,idxs]=sort(num_reactions);

end