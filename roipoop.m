
a0 = geometry.a0;

image_Ny = geometry.imageNy;
image_Nz = geometry.imageNz;

[Q, QQ, allK_q, allhkl] = generate_reduced_wavevectors(geometry);
%%
roic = [343 463 157 275];

roinz = roic(2)-roic(1)+1;
roiny = roic(4)-roic(3)+1;

tmp = reshape(QQ, [size(QQ,1) image_Ny image_Nz]);
roiQQ = tmp(:,roic(3):roic(4), roic(1):roic(2));
roiQQ = reshape(roiQQ, size(tmp,1), []);

tmp = reshape(Q', [size(Q,2) image_Ny image_Nz]);
roiQ = tmp(:,roic(3):roic(4), roic(1):roic(2));
roiQ = reshape(roiQ, size(tmp,1), [])';

%%

%peak_width = 0.005;

%mask = exp(-(roiQQ(1,:).^2+roiQQ(2,:).^2)/peak_width^2);

FWHM_kspace_bandwidth = 0.004;

sigma = FWHM_kspace_bandwidth/2/sqrt(2*log(2));

mask = exp(-(roiQQ(1,:).^2+roiQQ(2,:).^2)/2/sigma^2);

imagesc(reshape(mask, [roiny roinz]))


%%

Nz = 80;

r1z = 0:a0:a0*(Nz-1);
r2z = r1z + a0/2;


% atomic displacement from equilibrium in angstroms
u = -10*dxsamp*fliplr(cumsum(fliplr(strain(:,sampind)), 2));
%u = 10*dxsamp*cumsum(strain(:,sampind), 2);

rz0 = 0:a0/2:(Nz-1/2)*a0;
atom1index = 1:2:length(rz0);
atom2index = atom1index + 1;

r1z0 = rz0(atom1index);
r2z0 = rz0(atom2index);

[simdepth, simtmesh] = meshgrid(xsamp-xsamp(1), t(1:end-1));
[atomic_depth, atomictmesh] = meshgrid(rz0/10, t(1:end-1));

atomic_u = interp2(simdepth, simtmesh, u, atomic_depth, atomictmesh);

rz = bsxfun(@plus, atomic_u, rz0);
rz = bsxfun(@minus, rz, rz(:,1));
r1z = rz(:,atom1index);
r2z = rz(:,atom2index);

phi = squeeze(2*pi*roiQ(:,1))*a0/2 + squeeze(2*pi*roiQ(:,2))*a0/2;

Qz = squeeze(roiQ(:,3));

%%
for j = 1:size(rz, 1)
plot(rz0)
hold on
plot(rz(j,:), 'r')
hold off
drawnow
pause(0.01)
end

%%

ymax = 3;

for j = 1:size(rz, 1)
plot(xsamp-xsamp(1), u(j,:))
hold on
plot(rz0/10, atomic_u(j,:) , 'r')
hold off
ylim([-ymax ymax])
drawnow
pause(0.01)
end

%%

atl = 300;

%ef = sum(exp(1i*2*pi*Qz*(r1z0)) + exp(1i*2*pi*bsxfun(@plus, Qz*(r2z0), phi)), 2);

ef = sum( bsxfun(@times, exp(1i*2*pi*Qz*(r1z0)), exp(-r1z0/atl))...
    + bsxfun(@times, exp(1i*2*pi*bsxfun(@plus, Qz*(r2z0), phi)),exp(-r2z0/atl)) , 2);

Izoff = abs(ef).^2;

Ioff = reshape(Izoff.*mask', [roiny roinz]);

%imagesc(Iold)


%%
t0ind=45;
tendind = 401;

tsim = t(1:2:tendind);

Iz = zeros(length(tsim), length(Qz));


tic

Iz(1:(t0ind+1)/2,:) = bsxfun(@plus, Iz(1:(t0ind+1)/2,:), Izoff');



ctr = (t0ind+1)/2;

for tind = t0ind+2:2:tendind
    
    ctr = ctr+1;
    
    %ef = sum(exp(1i*2*pi*Qz*r1z(tind,:)) +...
    %    exp(1i*2*pi*bsxfun(@plus, Qz*r2z(tind,:), phi)), 2);
    
    ef = sum( bsxfun(@times, exp(1i*2*pi*Qz*(r1z(tind,:))), exp(-r1z(tind,:)/atl))...
    + bsxfun(@times, exp(1i*2*pi*bsxfun(@plus, Qz*(r2z(tind,:)), phi)),exp(-r2z(tind,:)/atl)) , 2);

    
    Iz(ctr,:) = abs(ef).^2;

end
toc

%%

beamwidth = 10;

tmp = -20:20;

beamgauss = exp(-tmp.^2/beamwidth^2);

%testcut = conv(I(:,50), beamgauss);

smearedIoff = conv2(Ioff, beamgauss', 'same');

imagesc(smearedIoff)

%plot(testcut)
%%
record_flag = 0; % set to 0 for NOT recording a video and 1 for recording

if record_flag == 1
    writerObj = VideoWriter('xray_simulation_smeared_run170.avi');
    writerObj.FrameRate = 10;
    writerObj.Quality = 100;
    open(writerObj);
    
end

I =  bsxfun(@times, Iz, mask);
I = reshape(I, [size(I,1) roiny roinz]);
smearedI = zeros(size(I));
%%
for j = 1:length(tsim)
    
    
    %smearedI(j,:,:) = conv2(squeeze(I(j,:,:)), beamgauss', 'same');
    %
    imagesc(squeeze(smearedI(j,:,:))-smearedIoff)
    
    axis image
    caxis([0, 45])
    
    title(['t=' num2str(tsim(j)) ' ps'])
    %colorbar
    drawnow
    pause(0.02)
    
    if record_flag == 1
        frame = getframe;
        writeVideo(writerObj,frame);
    end
    %}
    
end

if record_flag == 1
    close(writerObj);
end

%%


tmp = reshape(bsxfun(@times, Iz, mask), [size(Iz, 1) roiny roinz]);
sim_roi_int = sum(sum(tmp(:,37:65,36:82), 3), 2);
sim_roi_int = sim_roi_int - mean(sim_roi_int);
sim_roi_int = sim_roi_int/max(sim_roi_int);
plot(tsim, sim_roi_int)
xlabel('time delay (ps)')
ylabel('normalized roi intensity')
%plot(t(1:end-1), sum(sum(tmp(:,128:244,393:461), 3), 2))
hold all


%%

sim_roi_int = sum(bsxfun(@times, Iz, mask), 2);
%sim_roi_int = sum(sum(tmp(:,199:235,418:441), 3), 2);
sim_roi_int = sim_roi_int - mean(sim_roi_int);
sim_roi_int = sim_roi_int/max(sim_roi_int);
plot(tsim, sim_roi_int)
xlabel('time delay (ps)')
hold all

%%

maxint = max(max(smearedIoff));
maxroi = smearedIoff > 0.5*maxint;

roi_off = sum(smearedIoff(:).*maxroi(:));
roi_int = bsxfun(@times, reshape(smearedI, size(smearedI,1),[]), maxroi(:)');
roi_int = sum(roi_int,2);
roi_pc = 1e2*(roi_int-roi_off)./roi_off;
plot(tsim, roi_pc)
xlabel('time delay (ps)')
ylabel('% change')
hold all

%%

maxint = max(max(smearedIoff));
maxroi = smearedIoff > 0.5*maxint;
roi_off = smearedIoff(maxroi);

Iflat = reshape(smearedI,size(smearedI,1),[]);
roi = Iflat(:,maxroi(:));
roi_diff = bsxfun(@minus, roi, roi_off');
roi_pc = 1e2*bsxfun(@rdivide, roi_diff, roi_off');
plot(tsim, mean(roi_pc,2))
xlabel('time delay (ps)')
ylabel('% change')
hold all

%%

atl = 200;

%ef = sum(exp(1i*2*pi*Qz*(r1z0)) + exp(1i*2*pi*bsxfun(@plus, Qz*(r2z0), phi)), 2);

ef = sum( bsxfun(@times, exp(1i*2*pi*Qz*(r1z0)), exp(-r1z0/atl))...
    + bsxfun(@times, exp(1i*2*pi*bsxfun(@plus, Qz*(r2z0), phi)),exp(-r2z0/atl)) , 2);

Izoff = abs(ef).^2;

Ioff = reshape(Izoff.*mask', [image_Ny image_Nz]);

figure();
imagesc(Ioff)
axis(ax)
%%
Nz = 80;

r1z = 0:a0:a0*(Nz-1);
r2z = r1z + a0/2;

phi = squeeze(2*pi*Q(:,1))*a0/2 + squeeze(2*pi*Q(:,2))*a0/2;

Qz = squeeze(Q(:,3));

ef = sum(exp(1i*2*pi*Qz*r1z) + exp(1i*2*pi*bsxfun(@plus, Qz*r2z, phi)), 2);

Iz = abs(ef).^2;

I = reshape(Iz.*mask', [image_Ny image_Nz]);

imagesc(I)


%%

Nz = 40;

r1z = 0:a0/2:(Nz-1/2)*a0;
r2z = r1z + a0/2;

Qz = squeeze(Q(:,3));

ef = sum(exp(1i*2*pi*Qz*r1z), 2);
It = abs(ef).^2;

It = reshape(It, [512 512]);
imagesc(It)

%%

at_plane = 300;

ef = zeros(length(Qz), length(t)-1);
tic
for at_plane = 1:length(r1z0)
    
    ef = ef + exp(1i*2*pi*Qz*r1z(:,at_plane)') +...
        exp(1i*2*pi*bsxfun(@plus, Qz*r2z(:,at_plane)', phi));
    

end
toc