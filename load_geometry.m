function geometry = load_geometry(detector, material)

%% default parameters

beam_center        = [16 2231];

    alpha           = 0.4;
    chi             = -0;
    phi0            = -33.5 - 0.6;  % 22.8 ;% 68.764;
    A               = 1;
    B               = 0;
    CC               = 0;
    T               = 300;

image_num_horz_pixels = 512;
image_num_vert_pixels = 512;

%% Load the Material Parameters

% generate the material filename
material_fn=strcat('./gmtmatfiles/', material, '.mat');

% load the material structure
m=load(material_fn);


%% sample orientation

sn=[0;0;1]; % miller indicies of sample normal

v=m.convrealvecs*sn;

axis=cross(v, [0;0;1])/norm(v);

angle=asind(norm(axis));

if angle==0
    sorot=eye(3); %if the angle is zero, the rotation matrix is just identity
else
    sorot=rotationmat3D(angle, axis); % sample orientation rotation matrix
end


%% makes the Rotation matrix
Rot 	= eye(3);		% eye(3) is the identity

Rot		= rotationmat3D(phi0,[0 0 1])*Rot;		% rotation along sample normal [deg]
Rot		= rotationmat3D(alpha,[0 -1 0])*Rot;	% alpha = incidence angle [deg]
Rot		= rotationmat3D(chi,[1 0 0])*Rot;	% chi = rotation angle about x ray direction [deg]



%% Load detector parameters

detector_fn=strcat('./gmtmatfiles/', detector, '.mat');

detector_struct=load(detector_fn);

detector_struct.det_dist=195;

%% Construct Geometry

geometry = struct(  'a0',m.a0,...
                    'b0',m.b0,...
                    'c0',m.c0,...
                    'Nat',m.Nat,...
                    'Dmatrix_func',m.Dmatrix_func,...
                    'realvecs',m.realvecs,...
                    'convrealvecs',m.convrealvecs,...
                    'sorot',sorot,...
                    'primvects',m.primvects,...
                    'basis',m.basis,...
                    'Ws',m.Ws,...     % Debye-Waller 
                    'fs',m.fs,...     % Atomic factors
                    'Ms',m.Ms,...   % Ionic Masses in 10^-24 kg
                    'Rot',Rot,...
                    'lambda0',12.398/7,...
                    'imageNy',image_num_horz_pixels,...
                    'imageNz',image_num_vert_pixels,...
                    'detector',detector_struct,...
                    'beam_center',beam_center,...
                    'alpha',alpha,...
                    'chi',chi,...
                    'phi',phi0,...
                    'tdsA',A,...
                    'tdsB',B,...
                    'tdsC',CC,...
                    'temp',T);
                
                