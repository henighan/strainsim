function [strainstuff, expdet] = strainfunc_vthick(expdet, Te, Ti, fthick)


F0 = expdet.F0;
samp = expdet.samp;
sub = expdet.sub;

newendind = 1;

%rc = 0;
%tc = 1;
Temp0 = 300;

%F0 = 40; %laser fluence in mJ/cm^2
elpr0 = 1e-5; %electron pressure constant


% Thermal expansion coefficient and bulk modulus. These should really be
% added to the mat files.
sampbeta = 1.2e-5;
B = 1.7e11;

% Pulse arrives at t=0, dt is the spacing between time points in the strain
% simulation (temperature simulation may need smaller time steps, but this
% is done in tempsim). Times are in picoseconds. Note that dt determines
% dx. We want 1 time step to be 1 xstep. 
tstart = expdet.tstart;
tend = expdet.tend;
dt = expdet.dt;
t = expdet.t;

% xfront is the thickness of the cap layer, xback is where the sample ends
% and the substrate begins, xend is how far into the substrate we will
% simulate. 

%xfront = expdet.xfront;
%xback = expdet.xback;
%xend = expdet.xend;

% dt determines dx, 1 step in time is 1 step in x.
dxsub = expdet.dxsub;
dxsamp = expdet.dxsamp;

% make the x arrays for the cap layer, sample, and substrate
xcap = expdet.xcap;
xsamp = expdet.xsamp;
xsub = expdet.xsub;

lsampind = round(fthick/dxsamp);
xsamp = xsamp(1:lsampind);
xsub = xsub - xsub(1) + xsamp(end) + dxsamp;



Te = Te(:,1:lsampind);
Ti = Ti(:,1:lsampind);

%x = expdet.x;
x = [xcap xsamp xsub];

expdet.xsamp = xsamp;
expdet.xsub = xsub;
expdet.xcap = xcap;
expdet.x = x;

% reflection and transmission coefficients, calculated from the impedences
rc = (sub.imp - samp.imp)/(sub.imp + samp.imp);
tc = 2*sqrt(samp.imp*sub.imp)/(samp.imp + sub.imp);

%rc = -0.5;
%tc = sqrt(1-rc^2);

%newendind = 1;
%capind = 1:length(xcap);
capind = newendind:length(xcap);
sampind = length(xcap)+1:length(xcap) + length(xsamp);
subind = sampind(end)+1:length(x);




dTemp = [zeros(length(t), length(xcap)) Te(:,1:lsampind)-Temp0 zeros(length(t), length(xsub))];

Tewhole = dTemp+Temp0;

Tesq = Tewhole.^2;

dTesq = Tesq - Temp0^2;

dTesqdt = (dTesq(2:length(t),:) - dTesq(1:(length(t)-1),:))/dt;



%estr = ones(size(interptmesh'));
%estr = 4*Temp0*samp.CepT/3/(samp.rho*1e3*(1e3*samp.v)^2)*estr;
%estr = Tewholeinterp*samp.CepT/3/(samp.rho*1e3*(1e3*samp.v)^2);

%straini = 3*sampbeta*B/2/(samp.rho*1e3*(1e3*samp.v)^2)*...
    %(strainr + strainl + dTemp(2:end,:));

%% electron strain

%dTemp is Te-Temp0 (the change in electron temperature). Note that right 
%now the cap layer and substrate are not changing tempurature. 

%dTemp = [zeros(length(t), length(xcap)) Te-Temp0 zeros(length(t), length(xsub))];

% dTempdt is the rate of change of the ion temperature with respect to
% time. Remember time is in picoseconds here. 
%dTempdt = (dTemp(2:length(t),:) - dTemp(1:(length(t)-1),:))/dt;

sz = size(dTesqdt);

% Here we define the strain propogating to the right (strainr) and the
% strain propogating to the left (strainl).
strainr = zeros(sz);
strainl = zeros(sz);

%tic

% calculate the strain for the first timestep
strainr(1,:) = -0.5*dTesqdt(1,:)*dt;
strainl(1,:) = -0.5*dTesqdt(1,:)*dt;

% now loop through the other time steps and calculate the strain based on
% the change in temperature during that step and the propogating strain
% from the previous timestep.
for j = 2:sz(1)
    
        % right propogating strain from previous step in cap
        strainr(j,capind(2):capind(end)) = strainr(j-1,capind(1):capind(end-1));
        
        % right propogating strain from previous step in sample
        strainr(j,sampind(2):sampind(end)) = strainr(j-1,sampind(1):sampind(end-1));
        
        % right propogating strain from previous step in substrate
        strainr(j,subind(2):subind(end)) = strainr(j-1,subind(1):subind(end-1));
        
        % reflection off the free surface
        strainr(j,newendind) = strainr(j,newendind)-strainl(j-1,newendind);
        
        % reflection off cap back into sample
        %strainr(j,sampind(1)) = strainr(j,sampind(1)) + rc*strainl(j-1, sampind(1));
        strainr(j,sampind(1)) = strainr(j,sampind(1)) - strainl(j-1, sampind(1));
        %strainr(j,sampind(1)) = strainr(j,sampind(1)) + strainl(j-1, sampind(1));
        %strainr(j,sampind(1)) = strainr(j,sampind(1)) - rc*strainl(j-1, sampind(1));
        
        % transmission from cap into sample
        strainr(j,sampind(1)) = strainr(j,sampind(1)) + tc*strainr(j-1,capind(end));
        
        % transmission from sample into substrate
        strainr(j,subind(1)) = strainr(j,subind(1)) + tc*strainr(j-1,sampind(end));
        
        % right going strain from the temperature change in this step
        strainr(j,:) = strainr(j,:)-0.5*dTesqdt(j,:)*dt;
        
        % left propogating strain from previous step in cap
        strainl(j,capind(1):capind(end-1)) = strainl(j-1,capind(2):capind(end));
        
        % right propogating strain from previous step in sample
        strainl(j,sampind(1):sampind(end-1)) = strainl(j-1,sampind(2):sampind(end));
        
        % right propogating strain from previous step in substrate
        strainl(j,subind(1):subind(end-1)) = strainl(j-1,subind(2):subind(end));
        
        % reflection off sample back into cap
        strainl(j,capind(end)) = strainl(j,capind(end)) - rc*strainr(j-1,capind(end));
        
        % reflection off substrate back into sample
        strainl(j,sampind(end)) = strainl(j,sampind(end)) + rc*strainr(j-1,sampind(end));
        
        % transmission from sample into cap
        %strainl(j,capind(end)) = strainl(j,capind(end)) + tc*strainl(j-1, sampind(1));
    
        % left going strain from temperature change in this step
        strainl(j,:) = strainl(j,:)-0.5*dTesqdt(j,:)*dt;
        
end


estr = samp.CepT/3/(samp.rho*1e3*(1e3*samp.v)^2);
straine = estr*(strainr + strainl + dTesq(2:end,:));

%dTemp is Ti-Temp0 (the change in ionic temperature). Note that right 
%now the cap layer and substrate are not changing tempurature. 

dTemp = [zeros(length(t), length(xcap)) Ti-Temp0 zeros(length(t), length(xsub))];

% dTempdt is the rate of change of the ion temperature with respect to
% time. Remember time is in picoseconds here. 
dTempdt = (dTemp(2:length(t),:) - dTemp(1:(length(t)-1),:))/dt;

sz = size(dTempdt);

% Here we define the strain propogating to the right (strainr) and the
% strain propogating to the left (strainl).
strainr = zeros(sz);
strainl = zeros(sz);

%tic

% calculate the strain for the first timestep
strainr(1,:) = -0.5*dTempdt(1,:)*dt;
strainl(1,:) = -0.5*dTempdt(1,:)*dt;

% now loop through the other time steps and calculate the strain based on
% the change in temperature during that step and the propogating strain
% from the previous timestep.
for j = 2:sz(1)
    
        % right propogating strain from previous step in cap
        strainr(j,capind(2):capind(end)) = strainr(j-1,capind(1):capind(end-1));
        
        % right propogating strain from previous step in sample
        strainr(j,sampind(2):sampind(end)) = strainr(j-1,sampind(1):sampind(end-1));
        
        % right propogating strain from previous step in substrate
        strainr(j,subind(2):subind(end)) = strainr(j-1,subind(1):subind(end-1));
        
        % reflection off the free surface
        strainr(j,newendind) = strainr(j,newendind)-strainl(j-1,newendind);
        
        % reflection off cap back into sample
        %strainr(j,sampind(1)) = strainr(j,sampind(1)) + rc*strainl(j-1, sampind(1));
        strainr(j,sampind(1)) = strainr(j,sampind(1)) - strainl(j-1, sampind(1));
        %strainr(j,sampind(1)) = strainr(j,sampind(1)) + strainl(j-1, sampind(1));
        %strainr(j,sampind(1)) = strainr(j,sampind(1)) - rc*strainl(j-1, sampind(1));
        
        % transmission from cap into sample
        strainr(j,sampind(1)) = strainr(j,sampind(1)) + tc*strainr(j-1,capind(end));
        
        % transmission from sample into substrate
        strainr(j,subind(1)) = strainr(j,subind(1)) + tc*strainr(j-1,sampind(end));
        
        % right going strain from the temperature change in this step
        strainr(j,:) = strainr(j,:)-0.5*dTempdt(j,:)*dt;
        
        % left propogating strain from previous step in cap
        strainl(j,capind(1):capind(end-1)) = strainl(j-1,capind(2):capind(end));
        
        % right propogating strain from previous step in sample
        strainl(j,sampind(1):sampind(end-1)) = strainl(j-1,sampind(2):sampind(end));
        
        % right propogating strain from previous step in substrate
        strainl(j,subind(1):subind(end-1)) = strainl(j-1,subind(2):subind(end));
        
        % reflection off sample back into cap
        strainl(j,capind(end)) = strainl(j,capind(end)) - rc*strainr(j-1,capind(end));
        
        % reflection off substrate back into sample
        strainl(j,sampind(end)) = strainl(j,sampind(end)) + rc*strainr(j-1,sampind(end));
        
        % transmission from sample into cap
        %strainl(j,capind(end)) = strainl(j,capind(end)) + tc*strainl(j-1, sampind(1));

        % left going strain from temperature change in this step
        strainl(j,:) = strainl(j,:)-0.5*dTempdt(j,:)*dt;
        
end

% Thermal expansion coefficient and bulk modulus. These should really be
% added to the mat files.
% this is the linear coefficient of thermal expansion!!!!!
sampbeta = 1.2e-5;
B = 1.7e11;

% Now give our strain the right units (dimensionless) so it represents a
% fractional change in lattice spacing.
%straini = 3*sampbeta*B/(samp.rho*1e3*(1e3*samp.v)^2)*...
%    (strainr + strainl + dTemp(1:end-1,:));

straini = (3)*sampbeta*B/(samp.rho*1e3*(1e3*samp.v)^2)*...
    (strainr + strainl + dTemp(2:end,:));

staticstraini0 =...
    (3)*sampbeta*B/(samp.rho*1e3*(1e3*samp.v)^2)*dTemp(2:end,:)/F0;

staticstraine0 = estr*dTesq(2:end,:)/F0;


%strain = 3*sampbeta*B/(samp.rho*1e3*(1e3*samp.v)^2)*...
%    (strainr + strainl);

%straini = 3*sampbeta*B/(samp.rho*1e3*(1e3*samp.v)^2)*...
%    (dTemp(2:end,:));

%toc

strain = straini  + straine;
%strain = straini;
%strain = straine;

penprofile = exp(-(xsamp - xsamp(1))/expdet.pendepth);
penprofile = penprofile/sum(penprofile);
Timean = sum(bsxfun(@times, Ti, penprofile), 2);
Timean0 = (Timean - Temp0)/F0;



strainstuff = struct(  'strain',strain,...
                    'straini',straini,...
                    'straini0',straini/F0,...
                    'straine',straine,...
                    'straine0',straine/F0,...
                    'staticstraini0',staticstraini0,...
                    'staticstraine0',staticstraine0,...
                    'sampind',sampind,...
                    'xsamp',xsamp,...
                    'Timean',Timean,...
                    'Timean0',Timean0,...
                    'Temp0',Temp0,...
                    'pulsewidth',expdet.pulsewidth,...
                    'F0',F0,...
                    't',t);


end

