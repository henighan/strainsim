figure
hold all
ctr = 0;
clear lgd
for runnum = [165 187]
    ctr = ctr + 1;
    load(['roiSimMatfiles/run', num2str(runnum), 'roidata.mat'])
    plot(roi_delays, roi_pc)
    lgd{ctr} = ['run' num2str(runnum)];
    
end

legend(lgd)