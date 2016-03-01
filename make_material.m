function make_material(name)

m.CepT = 702; % electron heat capacity per Temperature in J m^-3 K^-2
m.Ci = 3.33e6; % ion heat capacity in J m^-3 K^-1
m.D = -1; % lattice diffusion coefficient in m^2/s
m.De = -1; % electron diffusion coefficient in m^2/s
m.kappa0 = 80; % room temperature thermal conductivity in W m^-1 K^-1
m.v = 5.13; % speed of sound in km/s (same as nm/ps)
m.rho = 7.87; % density in g cm^-3
m.imp = m.rho*m.v; % acoustic impedance
m.g = 41.9e17; % electron-phonon coupling constant in W m^-3 K^-1

save([name '.mat'], 'm');

end