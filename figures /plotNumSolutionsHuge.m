function plotNumSolutionsHuge()

%input: hugematrix EFMs and kOPt results structs 
load('hugematrixrunning_time.mat');
timesEFMs=running_times; timesEFMs(timesEFMs==0)=[]; clear runnning_times;
%load('/home/milievsk/Documents/copenhagen /MATLAB/pathConstructor/BuildMPS/BuildMPS/distance4Results/Succinate9case /hugematrixResults.mat')
load('hugematrixResults.mat');
timesKbest=running_times; clear running_times;

% time limit is 12 hours for EFMs
%numSolutionsEFMs=linspace(1,size(timesEFMs,1),size(timesEFMs,1));
%numSolutionsKbest=linspace(1,size(timesKbest,1),size(timesKbest,1));

%fig=figure;
%cummulative running times
cum_Kbest=cumsum(timesKbest/60);
cum_EFMs=cumsum(timesEFMs/60);

%plot(cum_EFMs);
%log scale 
semilogy(cum_EFMs);
hold on;
%stem(cum_EFMs,numSolutionsEFMs);
%plot(cum_Kbest);
%log scale
semilogy(cum_Kbest);
legend('EFMs','kOpt')
ylabel('Computation time (minutes)')
xlabel('Number of solutions')

view(-90,90)
set(gca,'ydir','reverse');
%print in tif format 
%print(fig,'NumSolutionsTIF','-dtiffn');
end

