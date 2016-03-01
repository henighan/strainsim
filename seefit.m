%function seefit(strainstuff, geometry, paramvecs, t0list, Ti0)

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

%i = 1;
ctr = 0;
for i= [1:3 5]
    ctr = ctr + 1;
    figure(i)
    load(['roiSimMatfiles/run', num2str(runs(i)), 'roidata.mat'])
    plot(roi_delays, roi_pc)
    grid on
    hold all
    load(['roiSimMatfiles/run', num2str(runs(i)), 'roisim.mat'])
    %plot(roi_delays(1:2:length(roi_delays)), roi_pc)
    plot(roi_delays + t0list(ctr), roi_pc)
    
end



%end
