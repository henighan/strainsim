
function [t0final, chi2final, npts] = optimizet0(geometry, strainstuff, datainfo, atl, t0start, ksbw, Ti)

%{
[~, raw_sim_pc, raw_sim_delays] = roisim4(geometry, strainstuff, atl);
disp([num2str(sum(isnan(raw_sim_pc)))...
    ' nans in raw_sim_pc at beginning of opt0'])
%}
%t0start = 0;
t0step = 0.001;
%{
t = strainstuff.t;
xsamp = strainstuff.xsamp;
sampind = strainstuff.sampind;
strain = strainstuff.strain;
%}

tstart = -0.5;
tend = 6;



data_delays = datainfo.data_delays;
data_pc = datainfo.data_pc;

sigma2 = std(data_pc(data_delays<tstart))^2;
data_pc = data_pc - mean(data_pc(data_delays<tstart));

data_pc = data_pc(((tend > data_delays) & (data_delays > tstart)));
data_delays = data_delays(((tend > data_delays) & (data_delays > tstart)));



%disp(['sigma2 is ' num2str(sigma2)])

%disp(['N is ' num2str(length(data_pc)-6)])

%{
disp(['length of data_delays is ' num2str(length(data_delays))])
disp(['length of data_pc is ' num2str(length(data_pc))])
disp(['number of nans in data_pc: ' num2str(sum(isnan(data_pc)))])
disp(['data_delays go from ' num2str(min(data_delays))...
    ' to ' num2str(max(data_delays))]);

figure(2)
plot(data_delays, data_pc)

disp(['atl is ' num2str(atl)])
%}
%[~, raw_sim_pc, raw_sim_delays] = roisim3(geometry, strainstuff, atl);
%[~, raw_sim_pc, raw_sim_delays] = roisim4(geometry, strainstuff, atl);
[~, raw_sim_pc, raw_sim_delays] = roisim_par_dwf(geometry, strainstuff, Ti, atl, ksbw);

%{
disp(['number of nans in raw_sim_pc: ' num2str(sum(isnan(raw_sim_pc)))])

disp(['raw_sim_delays go from ' num2str(min(raw_sim_delays))...
     ' to ' num2str(max(raw_sim_delays))])
%}

%sim_delays = raw_sim_delays(1:2:length(raw_sim_delays))+t0start;
sim_delays = raw_sim_delays+t0start;
sim_pc = interp1(sim_delays, raw_sim_pc, data_delays);

%{
disp(['length of sim_delays is ' num2str(length(sim_delays))])
disp(['length of sim_pc is ' num2str(length(sim_pc))])
disp(['number of nans in sim_pc is: ' num2str(sum(isnan(sim_pc)))])
%}
%%
chi2last = sum((sim_pc - data_pc').^2)/sigma2;
%chi2last = sum((sim_pc - data_pc').^2)/sigma2/(length(data_pc)-6);

%disp(['chi2last is ' num2str(chi2last)])
t0offset = t0start + t0step;

%sim_delays = raw_sim_delays(1:2:length(raw_sim_delays))+t0offset;
sim_delays = raw_sim_delays+t0offset;
sim_pc = interp1(sim_delays, raw_sim_pc, data_delays);
%{
disp(['length of sim_delays is ' num2str(length(sim_delays))])
disp(['length of sim_pc is ' num2str(length(sim_pc))])
%}

chi2 = sum((sim_pc - data_pc').^2)/sigma2;
%chi2 = sum((sim_pc - data_pc').^2)/sigma2/(length(data_pc)-6);
%disp(['chi2 is ' num2str(chi2)])

if chi2 < chi2last
    chi2last = chi2;
    t0offset = t0offset + t0step;
else
    t0step = -t0step;
    t0offset = t0offset + 2*t0step;
end

%sim_delays = raw_sim_delays(1:2:length(raw_sim_delays))+t0offset;
sim_delays = raw_sim_delays+t0offset;
sim_pc = interp1(sim_delays, raw_sim_pc, data_delays);

%disp(['length of sim_delays is ' num2str(length(sim_delays))])
%disp(['length of sim_pc is ' num2str(length(sim_pc))])

chi2 = sum((sim_pc - data_pc').^2)/sigma2;
%chi2 = sum((sim_pc - data_pc').^2)/sigma2/(length(data_pc)-6);
%disp(['chi2 is ' num2str(chi2)])

while chi2 < chi2last
    chi2last = chi2;
    
    t0offset = t0offset + t0step;
    
    %sim_delays = raw_sim_delays(1:2:length(raw_sim_delays))+t0offset;
    sim_delays = raw_sim_delays+t0offset;
    sim_pc = interp1(sim_delays, raw_sim_pc, data_delays);
    
    chi2 = sum((sim_pc - data_pc').^2)/sigma2;
    %chi2 = sum((sim_pc - data_pc').^2)/sigma2/(length(data_pc)-6);
%    disp(['chi2 is ' num2str(chi2)])
end


npts = length(data_pc);
t0final = t0offset - t0step;
disp(['t0final is ' num2str(t0final)]);
chi2final = chi2last;

%
%sum(isnan(sim_pc))
%sum(sim_pc)
figure(1)
plot(data_delays, sim_pc)
hold all
plot(data_delays, data_pc)
hold off
xlabel('delays (ps)')
ylabel('% change in diffracted intensity')
drawnow
%}
end
%{
plot(data_delays, data_pc)
sim_delays = raw_sim_delays(1:2:length(raw_sim_delays))+t0final;
sim_pc = interp1(sim_delays, raw_sim_pc, data_delays);
hold all
plot(data_delays, sim_pc)
%}