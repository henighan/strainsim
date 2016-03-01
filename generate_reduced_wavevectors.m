function [Q, QQ, allK_q, allhkl] = generate_reduced_wavevectors(geometry)
% This is the main script to set up the geometry parameters for the analysis 
% and calculation of TDS. It should be edited according to the experimental
% parameters. Never assume that the values used here are the right ones!
%
% This program should define the following parameters:
% 
% b1, b2, b3 		(global)
% these are the vector of the primitive unit cell in rec. space. You can
% use the definitions below (uncomment the parts you want) or define these
% on your own somewhere else, but make sure it is the PRIMITIVE basis.
% 
% a0				(may be global)
% lattice parameter (in Bi it's the length of a1, in ZB is the side of the cube)
% 
% Q					
% scattering vector for each pixel. make sure you edit 'det_kspace_proj.m'
% and set the detector size and number of pixels.
%
% detector_dist [mm]
%
% lambda0			(x-ray wavelength [Angstroms])
%
% Rot				(rotation of the lattice in experiment)
% 
% allhkl			(3 x size(Q,2) x size(Q,3) hkl for each pixel)
% this is used by functions like plotdata.m 
%
% allK_q			(same size as 'allhkl', rec.lat. vector for each pixel)
%
% QQ				(same size as 'allhkl', reduced wavevector for each pixel)

a0      = geometry.a0;
b0      = geometry.a0;
c0      = geometry.a0;

detector_dist= geometry.detector.det_dist;
x0      = geometry.beam_center;
alpha   = geometry.alpha;
chi     = geometry.chi;
phi0    = geometry.phi;
lambda0 = geometry.lambda0;

b1=geometry.primvects(:,1);
b2=geometry.primvects(:,2);
b3=geometry.primvects(:,3);

%a0      = 5.65325;		% GaAs lattice
%a0     = 5.658;		% Ge
%a0		= 6.479;		% InSb lattice constant
%a0		= 5.65325		% GaAs lattice
%a0		= 5.8687;		% InP lattice constant
%a0		= 5.97;			% SmS
%a0     = 5.658;		% Ge
%a0		= 5.431;		% Si

%a=4.785; c=12.991		% sapphire structure (hexagonal)

%a0		= 4.7236;		% Bi <------- Bismuth
%

%% useful to generate the BZ knowing b1, b2, b3
% first define the BZ

%{
b1= [-1 1 1]'/a0;
b2= [1 -1 1]'/a0;
b3= [1 1 -1]'/a0;

%%
%}
%lambda0 = 12.398/10;		% x-ray wavelength

Q = det_kspace_proj(geometry);

% NOTE: Rot should be the rotation matrix that rotates the crystal axes
% (a1,a2,a3) and it's inverse (Rot') rotates Q

Rot 	= geometry.Rot;

%disp('b1,b2,b3 axes in the rotated reference frame:')
%disp([(Rot*b1), (Rot*b2),(Rot*b3)])
%}

%% bismuth lattice and BZ
%{
ang=asin(2/sqrt(3)*sin(57.35/180*pi/2));
a1 = [a0*sin(ang) , 0, a0*cos(ang)];
a2 = [-a0*sin(ang)/2, a0*sin(ang)*sqrt(3)/2, a0*cos(ang)];
a3 = [-a0*sin(ang)/2, -a0*sin(ang)*sqrt(3)/2, a0*cos(ang)];

b1 = cross(a2,a3)'./(a1*cross(a2,a3)');
b2 = cross(a3,a1)'./(a1*cross(a2,a3)');
b3 = cross(a1,a2)'./(a1*cross(a2,a3)');
%}

%{
b1 = Rot*b1';
b2 = Rot*b2';
b3 = Rot*b3';
%}

%%
count = 0;
verts = [];
% makes G point and the closest first neighbors. Using these vertices
% we then use voronoin() to construct the Wigner-Seitz cell of [0 0 0].
% i.e. the Brillouin zone.

% BZ contains the vertices of the Brillouin zone
%
NN=1;
for x= -NN:NN
	for y= -NN:NN
		for z= -NN:NN
			count = count +1;
			verts(:,count) = Rot*(x*b1+y*b2+z*b3);
		end
	end
end

%[BZconst, C] = voronoin(verts');		% finds the Voronoi cell of [0 0 0]
%BZ = BZconst(C{14},:);
%tri = delaunayn(BZ(C{14},:),{'Qt','Qbb','Qc','Qz'});

%}

%% loop over  a bunch of rec.lat. points close the the Ewald sphere \pm 1/a0
% roughly look for hkl points that are close to the Ewald sphere and may
% show up in the detector. we look in a band +- 1/a0 around the sphere.
% 
count = 0;
verts = [];
verts_hkls = [];

NN=3;
for x= -NN:NN
	for y= -NN:NN
		for z= -NN:NN
			Ghkl = Rot*(x*b1+y*b2+z*b3);
			% look for points near the sphere and positive hemisphere
			if (	(norm(Ghkl+1/lambda0*[1 0 0]') < 1/lambda0 + 1/a0) && ...
					(norm(Ghkl+1/lambda0*[1 0 0]') > 1/lambda0 - 1/a0))
				count = count+1;
				verts(:,count) = Ghkl;		% if it is in the region, keep it
				verts_hkls(:,count) = [x, y, z]'; % and keep the hkl indices
			end
		end
	end
end

%% for each Q point, searches which of verts(:) is closest then saves the
%% indices. Then ind(i) is the index of verts closest to Q(:,i)

Q=reshape(Q,[],3);
%[ind, D] = knnsearch(verts',Q);	% does a brute force linear search
[ind, D] = knnsearchFEX(Q,verts');	

%save the closest rec.lat. vector to each pixel
allK_q = verts(:,ind);
%find the reduced wavevector
QQ = Q(:,:)' - allK_q;					% Reduced wavevector [A^-1]

allhkl = verts_hkls(:,ind);				% also the hkl indices

%reshape to the image shape to use with other scripts
%QQ = reshape(QQ,[3 size(Q,2) size(Q,3)]);
%allhkl = reshape(allhkl,[3 size(Q,2) size(Q,3)]);
%allK_q = reshape(allK_q,[3 size(Q,2) size(Q,3)]);


%%
% this section calculates the projection of Q onto QQ (scattering vector vs
% the reduced wavevector) which is a measure of the sensitivity to L/T 
% (long. or transverse) polarization.
% 
%{
for j=1:Ny
	for i=1:Nz
		long_and_trans(i,j)= abs(Q(:,i,j)'*QQ(:,i,j)/norm(Q(:,i,j))/norm(QQ(:,i,j)));
	end
end
%}

%% 
% this was a sophisticated way to find nearest neighbors but not worth it
% (slower) than brute force search. you shoudn't use it.
%{

tnn = zeros([size(Q,2) size(Q,3)]);

for ii=1:size(verts,2)
	tn = tsearchn(BZ(C{14},:),tri,Q(:,:)' - repmat(verts(:,ii)',[size(Q,2)*size(Q,3) 1]));
	tn(~isnan(tn))=ii;
	tn(isnan(tn))=0;
	tnn = reshape(tn, [size(Q,2), size(Q,3)])+tnn;
%	if ()
%		imagesc(tnn);axis image; 
%		title(num2str(ii))
%	end
end

imagesc(tnn);axis image; 

%}
