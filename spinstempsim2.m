function [Ti_out, Te_out, Ts_out] = spinstempsim2(material, expdet, F)


%[Ti, Te, Ts] = spinstempsim2(samp, expdet.xsamp, expdet.t, F0);

input_x = expdet.xsamp;
input_t = expdet.t;

%A = pi*(10e-6)^2; % area in m^-2
%Q = 0.1*A; % pulse energy in joules
%F = 50; %laser fluence in J m^-2
%F = 7; %laser fluence in J m^-2
R = 0.6; % optical reflection coefficient at surface
CepT = material.CepT; % electron heat capacity per unit temperature in J m^-3 K^-2
Ci = material.Ci; % ion heat capacity in J m^-3 K^-1
%De = material.De; % electron diffusion coefficient in m^2/s

g = material.g; % Electron-phonon coupling constant in W m^-3 K^-1
%g = 1000*material.g; % Electron-phonon coupling constant in W m^-3 K^-1
kappa0 = material.kappa0;
%kappa0 = 0;
Temp0 = 300;

% I'm ballparking this based on the measurements of Carpene.
%ges = 3*material.g; % Electron-spin coupling constant in W m^-3 K^-1.
ges = 0; % Electron-spin coupling constant in W m^-3 K^-1.
CspT = 7.5; %spin heat capacity per T^3/2 in J m^-3 K^-5/2

% For the simluation to be stable, it seems we need much shorter time-steps
% than what is used for the strain. 
tstart = input_t(1);
tend = input_t(end);
dt = 0.00001;
t = tstart:dt:tend;

% laser pulsewidth and penetration depth
%pulsewidth = 0.05; % in ps
%pendepth = 16; % in nm
pulsewidth = expdet.pulsewidth;
pendepth = expdet.pendepth;


targetdx = 0.1;
xstart = input_x(1);
xend = input_x(end);
nxsteps = round((xend - xstart)/targetdx);
dx = (xend - xstart)/nxsteps;

x = xstart:dx:xend;

%dx = x(2) - x(1);

% note that the integral of pulse profile is equal to 1 for infinite
% tickness film
pulseprofile = exp(-(t/pulsewidth).^2)/sqrt(pi)/pulsewidth;

penetrationprofile = exp(-(x-x(1))/pendepth);


%%
Te = zeros(length(t), length(x));
Ti = zeros(length(t), length(x));
Ts = zeros(length(t), length(x));

Te(1,:) = Temp0;
Ti(1,:) = Temp0;
Ts(1,:) = Temp0;


tic
for j = 2:length(t)
    
    Ce = CepT*Te(j-1,:);
    Cs = CspT*(Ts(j-1,:).^(3/2));
    
    Teright = circshift(Te(j-1,:)', -1);
    Teright(end) = Teright(end-1);
    Teleft = circshift(Te(j-1,:)', 1);
    Teleft(1) = Teleft(2);
    
    Tiright = circshift(Ti(j-1,:)', -1);
    Tiright(end) = Tiright(end-1);
    Tileft = circshift(Ti(j-1,:)', 1);
    Tileft(1) = Tileft(2);
    
    d2Tedx2 = (Teright' + Teleft' - 2*Te(j-1,:))/dx^2;
    
    dTedx = (Teright' - Teleft')/(2*dx);
    dTedx(1) = 0;
    dTedx(end) = 0;
    
    dTidx = (Tiright' - Tileft')/(2*dx);
    dTidx(1) = 0;
    dTidx(end) = 0;
    
    Te(j,:) = Te(j-1,:) +...
        dt*(1-R)*F*pulseprofile(j)*penetrationprofile./(Ce*(pendepth*1e-9)) + ...
        (1e6)*dt*kappa0*d2Tedx2./Ce + (1e6)*dt*kappa0*(dTedx.^2)./(Ce.*Ti(j-1,:)) + ...
        -(1e6)*dt*kappa0*Te(j-1,:).*dTedx.*dTidx./(Ce.*Ti(j-1,:).^2) + ...
        -(1e-12)*dt*g*(Te(j-1,:) - Ti(j-1,:))./Ce + ...
        -(1e-12)*dt*ges*(Te(j-1,:) - Ts(j-1,:))./Ce;
    
    
    Ti(j,:) =  Ti(j-1,:) + (1e-12)*dt*g*(Te(j-1,:) - Ti(j-1,:))/Ci;
    Ts(j,:) =  Ts(j-1,:) + (1e-12)*dt*ges*(Te(j-1,:) - Ts(j-1,:))./Cs;
    
    Te(j,1) = mean(Te(j,1:2));
    Te(j,2) = Te(j,1);
    Te(j,end) = mean([Te(j,end) Te(j,end-1)]);
    Te(j,end-1) = Te(j,end);
    
end
toc

[xmsh, tmsh] = meshgrid(x, t);

[xoutmsh, toutmsh] = meshgrid(input_x, input_t);

Ti_out = interp2(xmsh, tmsh, Ti, xoutmsh, toutmsh);
Te_out = interp2(xmsh, tmsh, Te, xoutmsh, toutmsh);
Ts_out = interp2(xmsh, tmsh, Ts, xoutmsh, toutmsh);

end
