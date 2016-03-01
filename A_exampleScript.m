% This is a quick demo of what this code is meant to do - fit the
% diffraction from a thin film excited by a short laser pulse


%% Initialization

% here we load in a matlab struct which contains a number of experimental
% and simulation details. Things like the film thickness, range of times
% we'll simulate over, the timesteps, etc.
expdet = load_expdet();

% in this particular example, I'm fixing the electron-phonon coupling
% constant to be fixed to the value below (in W/m^3/K).
expdet.samp.g = 1e18;

% normally I would initiate parallel computing here to speed things up
%matlabpool('open', 4)

% define the starting guess for the laser fluence
F0 = expdet.F0;

%% Temperature simulation (Two Temperature Model)

% here we simulate the spatial and temporal evolutions of the electronic
% and lattice temperatures by numerically solving the nonlinear
% heat-transfer differential equations
%[Ti, Te, Ts] = spinstempsim2(expdet.samp, expdet, F0);

% then save these temperature profiles so we don't have to calculate them
% again in the future if we re-run the fit
%save(['~/strainsim2/fitoutput/tempsgammaefit.mat'], 'Te', 'Ti', 'Ts');


% instead of simulating, I'm going to load the already simulated and saved
% temperatures to save us some time
load('./fitoutput/tempsgammaefit.mat');

% plot a visualization of the electron and lattice temperatures
seetemp(expdet, Te, Ti)
%% Strain simulation

% here we load in a struct (which is called geometry for historical reasons
% in my research group's suite of code) which contains details relevant for
% calculating the x-ray diffraction patterns from the sample. These include
% information about the sample's crystal lattice and its orientation, the
% position of the area detector and how many pixels it has, etc.
geometry = load_geometry('CSPAD', 'Fe');

fthick = 22.6; % this is the film thickness

% simulate the strain caused by the temperatures simulated above given our
% starting parameters
[strainstuff, expdet] = strainfunc_vthick(expdet, Te, Ti, fthick);

strainstuff.fthick = fthick;
Ti0 = (Ti-300)/F0;
Ti0 = Ti0(:,1:length(expdet.xsamp));

% visualize the simulated strain in the film and the substrate
seestrain(expdet, strainstuff)
%% Calculate Diffraction and fit to measured diffraction

%elstrcoeff = 0;
atl = 660;
phioffset = -0.66;
elstrcoeff0 = 3.6;
%ksbw0 = 0.006;

% these will be the free parameters in the fit, the fluence of the laser
% pulse, our errors in the sample orientation angles, and the electron
% strain coefficient, which is a physical parameter that specifies the
% thermal expansion coefficient of the electrons
paramvec0 = [F0 phioffset atl elstrcoeff0];

chi2handle = @gammaeBiggerchi2;
paramvecout = paramvec0;
options = optimset('MaxFunEvals',1e10, 'MaxIter', 1e10);
%

% here we minimize the reduced chi^2 to obtain the best fit values of the
% above parameters. On each iteration, we adjust the input strain profile
% (simulated above) in correspondnce with changes in our free parameters,
% re-calculate the diffraction, and then compare this to the measured
% diffraction to calculate the reduced chi^2. Here we are minimizing the
% reduced chi^2.
%tic
paramvecout = ...
    fminsearch(@(paramvec) gammaeBiggerchi2(strainstuff, geometry, Ti0, paramvec), ...
    paramvec0, options)
%toc

%}

% return the reduced chi^2 and time0 offsets of our final fit
[chi2, t0list] = gammaeBiggerchi2(strainstuff, geometry, Ti0, paramvecout);


% save the fitted parameters
fitstuff = struct(  'paramvecout',paramvecout,...
                    'chi2',chi2,...
                    't0list',t0list);
                
save(['~/strainsim2/fitoutput/gammaefreeG1.mat'], '-struct', 'fitstuff');
%}
%end

%% Visualize Fitted profiles against measured profiles


load('fitoutput/gammaefreeG1.mat')

paramvecs = paramvecout;

F = paramvecs(1);
%elstrcoeff = paramvecs(2);
phioffset = paramvecs(2);
atl = paramvecs(3);
%ksbw = paramvecs(4);
elstrcoeff = paramvecs(4);
ksbw = 0.006;

Ti = strainstuff.Temp0 + F*Ti0;

%elstrcoeff = 3.6;
%elstrcoeff = 0;
strainstuff.strain = F*strainstuff.straini0 + F*elstrcoeff*strainstuff.straine0;
strainstuff.Timean = strainstuff.Temp0 + F*strainstuff.Timean0;

philist = [-32 -42 -33.5 -36 -30 -29];
runs = [180 184 168 166 175 183];



geometry.alpha = 0.4;
%i=3;
%

for i=[1:3 5]
    
    phi = philist(i) + phioffset;
    
    %Rotation matrix for this phi
            Rot 	= eye(3);		% eye(3) is the identity
            
            Rot		= rotationmat3D(phi,[0 0 1])*Rot;		% rotation along sample normal [deg]
            Rot		= rotationmat3D(geometry.alpha,[0 -1 0])*Rot;	% alpha = incidence angle [deg]
            Rot		= rotationmat3D(geometry.chi,[1 0 0])*Rot;	% chi = rotation angle about x ray direction [deg]
            
            geometry.Rot = Rot;
    
    [roi_off, roi_pc, roi_delays] = roisim_par_dwf(geometry, strainstuff, Ti, atl, ksbw);
    save(['roiSimMatfiles/run', num2str(runs(i)), 'roisim.mat'], 'roi_off', 'roi_pc', 'roi_delays');
    %roioffsimlist(i) = roi_off;
    
end
%}
%
lw = 2; % plot linewidth
fs = 20; % plot fontsize

%i = 1;
ctr = 0;
for i= [1:3 5]
    ctr = ctr + 1;
    figure(i)
    load(['roiSimMatfiles/run', num2str(runs(i)), 'roidata.mat'])
    plot(roi_delays, roi_pc, 'linewidth', lw)
    set(gca, 'fontsize', fs)
    grid on
    hold all
    load(['roiSimMatfiles/run', num2str(runs(i)), 'roisim.mat'])
    %plot(roi_delays(1:2:length(roi_delays)), roi_pc)
    plot(roi_delays + t0list(ctr), roi_pc, 'linewidth', lw)
    xlabel('Delays (ps)')
    ylabel('% Change in Diffracted Intensity')
    legend('Data', 'Fit')
    
end



