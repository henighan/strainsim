function [roi_off, roi_pc, roi_delays] = roisim_par_dwf(geometry, strainstuff, Ti, atl, ksbw)

% this function simulates the time-resolved diffraction from the film given
% the geometry and strain

% we'll use this to calculate the debye-waller attenuation
Timean = strainstuff.Timean;
%disp([num2str(sum(isnan(Timean(:)))) ' nans in Timean']);

fthick = strainstuff.fthick;

tind = 1:1:length(strainstuff.t);
%disp(['lengthy of tind is ' num2str(length(tind))])
%disp(['size of strain is ' num2str(size(strainstuff.strain))])

strain = strainstuff.strain(tind(1:end-1),:);
sampind = strainstuff.sampind;
t = strainstuff.t(tind);
xsamp = strainstuff.xsamp;

%{
disp([num2str(sum(isnan(strain(:)))) ' nans in strain'])
disp([num2str(sum(isnan(sampind(:)))) ' nans in sampind'])
disp([num2str(sum(isnan(t(:)))) ' nans in t'])
disp([num2str(sum(isnan(xsamp(:)))) ' nans in xsamp'])
%}

a0 = geometry.a0;

image_Ny = geometry.imageNy;
image_Nz = geometry.imageNz;

imdim = [image_Ny image_Nz];

% Find the scattering vector for each pixel
[Q, QQ, ~, ~] = generate_reduced_wavevectors(geometry);

% define the sample normal vector
samp_normal = rotationmat3D(geometry.alpha, [0 -1 0])*[0 0 1]';

% find projection of each pixels scattering vector onto the sample normal
QQdotSN = sum(bsxfun(@times, QQ, samp_normal), 1);

% find the vector connecting the scattering vector and the
% projection of the scattering vector onto the sample normal.
% This vector will be perpendicular to the sample normal.
vector_to_trunc_rod = QQ - kron(samp_normal, QQdotSN);

% Find the length of the above vector
dist_to_trunc_rod = sum(vector_to_trunc_rod.^2, 1);

% reshape this 1d list of vector lenghts into a 2d image
dist_to_trunc_rod_im = reshape(dist_to_trunc_rod, imdim);

% Find the indicies of pixel which is closest to the truncation rod.
minind = find(dist_to_trunc_rod_im==(min(min(dist_to_trunc_rod_im))));

[j, i] = ind2sub(imdim, minind);

% position in mm of bottom left corner (blc) of detector, where the beam is
% going through (y,z) = (0,0)
detblcyz = bctodetyz(geometry, geometry.beam_center) +...
    [geometry.detector.det_size_horz/2 -geometry.detector.det_size_vert/2];

% distance of truncation rod pixel from bottom left corner of detector, in
% mm
yblc = -i*geometry.detector.det_size_horz/image_Ny;
zblc = (image_Nz - j)*geometry.detector.det_size_vert/image_Nz;


% position in mm of truncation rod pixel. 
ymm = yblc + detblcyz(1);
zmm = zblc + detblcyz(2);

%distance of truncation rod pixel from beam in the yz plane
yzrmm = sqrt(ymm^2 + zmm^2);

% truncation rod scattering angle 
%(angle between incoming and scattered xrays)
scatangle = atand(yzrmm/geometry.detector.det_dist);



% we'll only calculate pixels in an roi around the truncation rod to save
% time
roisize = 40;

% the corners of the roi (region of interest)
roic = [i-roisize i+roisize j-roisize j+roisize];

%number of roi pixels in z and y directions
roinz = roic(2)-roic(1)+1;
roiny = roic(4)-roic(3)+1;

%reduced scattering vectors in our roi
tmp = reshape(QQ, [size(QQ,1) image_Ny image_Nz]);
roiQQ = tmp(:,roic(3):roic(4), roic(1):roic(2));
roiQQ = reshape(roiQQ, size(tmp,1), []);

% scattering vectors in our roi
tmp = reshape(Q', [size(Q,2) image_Ny image_Nz]);
roiQ = tmp(:,roic(3):roic(4), roic(1):roic(2));
roiQ = reshape(roiQ, size(tmp,1), [])';

% the bandwidth in k-space of the FEL beam (in inverse angstroms. FEL
% bandwidth ~50 eV)
%FWHM_kspace_bandwidth = 0.006;
FWHM_kspace_bandwidth = ksbw;
sigma = FWHM_kspace_bandwidth/2/sqrt(2*log(2));

% give our truncation rod some width in the x and y directions to account
% for the finite bandwidth of the FEL. (crystal mosaicity could also be a
% factor here). Mask out other pixels.
mask = exp(-(roiQQ(1,:).^2+roiQQ(2,:).^2)/2/sigma^2);

%imagesc(reshape(mask, [roiny roinz]))

% number of unit cells in the sample normal direction
%Nz = 80;
%Nz = 89;
Nz = round(10*fthick/geometry.a0);
%disp(['Nz is ' num2str(Nz)])

%r1z = 0:a0:a0*(Nz-1);
%r2z = r1z + a0/2;

dxsamp = xsamp(2)-xsamp(1);

% atomic displacement from equilibrium in angstroms
u = -10*dxsamp*fliplr(cumsum(fliplr(strain(:,sampind)), 2));

% equilibrium positions of atomic planes
rz0 = 0:a0/2:(Nz-1/2)*a0;

% here I'm defining two indicies, each representing one of the two atomic
% planes in the conventional unit cell
atom1index = 1:2:length(rz0);
atom2index = atom1index + 1;

% equilibrium positions of the two sets of atomic planes
r1z0 = rz0(atom1index);
r2z0 = rz0(atom2index);

% we'll keep track of the displacement of each atomic plane at each
% timestep.
[simdepth, simtmesh] = meshgrid(xsamp-xsamp(1), t(1:end-1));
[atomic_depth, atomictmesh] = meshgrid(rz0/10, t(1:end-1));

atomic_u = interp2(simdepth, simtmesh, u, atomic_depth, atomictmesh);
atomic_Ti = interp2(simdepth, simtmesh, Ti(1:end-1,:), atomic_depth, atomictmesh);

% z positions of each atomic plane
rz = bsxfun(@plus, atomic_u, rz0);
rz = bsxfun(@minus, rz, rz(:,1));
% z positions of the two atomic planes associated with each conventional
% unit cell along the sample normal
r1z = rz(:,atom1index);
r2z = rz(:,atom2index);
% ionic temperature for each atomic plane
Ti1z = atomic_Ti(:,atom1index);
Ti2z = atomic_Ti(:,atom2index);

% This phi does NOT pertain to the sample orientation
% its the phase difference between x-rays scattered from the two different
% planes in the conventional unit cell due to the relative atomic displacement
% perpendicular to the sample normal
phi = squeeze(2*pi*roiQ(:,1))*a0/2 + squeeze(2*pi*roiQ(:,2))*a0/2;

Qz = squeeze(roiQ(:,3));

% first we calculat the diffraction without laser pumping. ef is the
% calculated electric field, given by the equations in the supplemental
% information of the paper. 
ef = sum( bsxfun(@times, exp(1i*2*pi*Qz*(r1z0)), exp(-r1z0/atl))...
    + bsxfun(@times, exp(1i*2*pi*bsxfun(@plus, Qz*(r2z0), phi)),exp(-r2z0/atl)) , 2);

% Calculate the diffracted intensity given the electric field
Izoff = (debyewaller(scatangle/2, geometry.lambda0, Timean(1))*abs(ef)).^2;

Ioff = reshape(Izoff.*mask', [roiny roinz]);

inputdt = t(2)-t(1);

targetdt = 0.04;

nskip = round(targetdt/inputdt);

%t0ind=45;
[~,t0ind] = min((t+strainstuff.pulsewidth).^2);

t0ind = t0ind + 1 - mod(t0ind,nskip);

tendind = length(t)-nskip;

tsim = t(1:nskip:tendind);

Iz = zeros(length(tsim), length(Qz));

%%
tic

Iz((1:(t0ind+nskip-1)/nskip),:) = bsxfun(@plus, Iz((1:(t0ind+nskip-1)/nskip),:), Izoff');


ctr0 = (t0ind+nskip-1)/nskip;
tinds = t0ind+nskip:nskip:tendind;

% here we loop over each timestep and calculate the summed electric field
% and intensity.
parfor i = 1:length(tinds)
    
    %ctr = ctr0+i;
    
    %ef = sum(exp(1i*2*pi*Qz*r1z(tind,:)) +...
    %    exp(1i*2*pi*bsxfun(@plus, Qz*r2z(tind,:), phi)), 2);
    
    %ef = sum( bsxfun(@times, exp(1i*2*pi*Qz*(r1z(tinds(i),:))), exp(-r1z(tinds(i),:)/atl))...
    %+ bsxfun(@times, exp(1i*2*pi*bsxfun(@plus, Qz*(r2z(tinds(i),:)), phi)),exp(-r2z(tinds(i),:)/atl)) , 2);
    
    
    ef = sum( bsxfun(@times, exp(1i*2*pi*Qz*(r1z(tinds(i),:))), ...
        debyewaller(scatangle/2, geometry.lambda0, Ti1z(tinds(i),:)).*exp(-r1z(tinds(i),:)/atl))...
        + bsxfun(@times, exp(1i*2*pi*bsxfun(@plus, Qz*(r2z(tinds(i),:)), phi)), ...
        debyewaller(scatangle/2, geometry.lambda0, Ti2z(tinds(i),:)).*exp(-r2z(tinds(i),:)/atl)) , 2);
    
    
    
    %Iz(ctr0+i,:) = debyewaller(scatangle/2, geometry.lambda0, Timean(tinds(i)))*abs(ef).^2;
    
    Iz(ctr0+i,:) = abs(ef).^2;
    %test = debyewaller(scatangle/2, geometry.lambda0, Timean(tinds(i)))*abs(ef).^2;
end
toc

% Here we consider the finite size of the xray beam, which was about 10
% pixels. We convolve the diffracted x-ray image with a gaussian of this
% size to approximate the smearing effect from the finite sized beam

beamwidth = 10;
%beamwidth = 3;

tmp = -20:20;

beamgauss = exp(-tmp.^2/beamwidth^2);

%testcut = conv(I(:,50), beamgauss);

smearedIoff = conv2(Ioff, beamgauss', 'same');

%imagesc(smearedIoff)


I =  bsxfun(@times, Iz, mask);
I = reshape(I, [size(I,1) roiny roinz]);
smearedI = zeros(size(I));


for j = 1:length(tsim)
    
    
    smearedI(j,:,:) = conv2(squeeze(I(j,:,:)), beamgauss', 'same');
    
    
end

% Finally, we output the % change in the diffracted intensity as a function
% of time within a region of interest

maxint = max(max(smearedIoff));
maxroi = smearedIoff > 0.5*maxint;

roi_off = sum(smearedIoff(:).*maxroi(:));
roi_int = bsxfun(@times, reshape(smearedI, size(smearedI,1),[]), maxroi(:)');
roi_int = sum(roi_int,2);
roi_pc = 1e2*(roi_int-roi_off)./roi_off;

%roi_delays = t(1:end-1);
roi_delays = tsim;
roi_off = mean(smearedIoff(maxroi));
%{
plot(tsim, roi_pc)
xlabel('time delay (ps)')
ylabel('% change')
hold all
%}
end