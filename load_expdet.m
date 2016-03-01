
function expdet = load_expdet()

% these are binary files in which I stored relevant information
% about the film (Fe, iron) and substrate (Magnesium Oxide, MgO) materials

load('Fe2.mat')
samp = m;
load('MgO.mat')
sub = m;

F0 = 15; %laser fluence in J/m^2

% Pulse arrives at t=0, dt is the spacing between time points in the strain
% simulation (temperature simulation may need smaller time steps, but this
% is done in tempsim). Times are in picoseconds. Note that dt determines
% dx. We want 1 time step to be 1 xstep. 
tstart = -1;
tend = 8;
dt = 0.01;
t = tstart:dt:tend;

% here x is the spatial coordinate normal to the film
% t is the tmporal coordinate

% xfront is the thickness of the cap layer, xback is where the sample ends
% and the substrate begins, xend is how far into the substrate we will
% simulate. 
%{
xfront = 5.2;
xback = 30.9;
xend = 304;
%}
%{
xfront = 10;
xback = 33;
xend = 303;
%}
xfront = 4.5;
%xback = 27.5;
xback = 27.55;
xend = 297.5;


% dt determines dx, 1 step in time is 1 step in x.
dxsub = sub.v*dt;
dxsamp = samp.v*dt;

% make the x arrays for the cap layer, sample, and substrate
xcap = 0:dxsub:xfront;
xsamp = (xcap(end)+dxsub):dxsamp:xback;
xsub = (xsamp(end)+dxsamp):dxsub:xend;

xfront = xsamp(1);
xback = xsub(1);

% truncate off machine-sized errors. matlab.....
xcap = round((10^10)*xcap)/(10^10);
xsamp = round((10^10)*xsamp)/(10^10);
xsub = round((10^10)*xsub)/(10^10);

x = [xcap xsamp xsub];

%laser pulse length and penetration depth
pulsewidth = 0.05; %pulse length in picoseconds
pendepth = 17.46; %penetration depth in nm, calculated from the Fresnel
% equations for Fe

expdet = struct(  'tstart',tstart,...
                    'tend',tend,...
                    'dt',dt,...
                    't',t,...
                    'F0',F0,...
                    'pulsewidth',pulsewidth,...
                    'pendepth',pendepth,...
                    'xfront',xfront,...
                    'xback',xback,...
                    'xend',xend,...
                    'dxsub',dxsub,...
                    'dxsamp',dxsamp,...
                    'xcap',xcap,...
                    'xsamp',xsamp,...
                    'xsub',xsub,...
                    'x',x,...
                    'sub',sub,...
                    'samp',samp);
                
end
