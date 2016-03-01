function geometry = changephi(geometry, phi)

%phi = geometry.phi + step;
alpha = geometry.alpha;

Rot 	= eye(3);		% eye(3) is the identity

Rot		= rotationmat3D(phi,[0 0 1])*Rot;		% rotation along sample normal [deg]
Rot		= rotationmat3D(geometry.alpha,[0 -1 0])*Rot;	% alpha = incidence angle [deg]
Rot		= rotationmat3D(geometry.chi,[1 0 0])*Rot;	% chi = rotation angle about x ray direction [deg]

geometry.Rot = Rot;

geometry.phi = phi;

end