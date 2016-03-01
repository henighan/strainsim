
%load('spinlesstemps15');
load('simtempstouse');

geometry = load_geometry('CSPAD', 'Fe');

expdet = load_expdet();

F0 = 15;
fthick = 22.6;
[strainstuff, expdet] = strainfunc_vthick(expdet, Te, Ti, fthick);
strainstuff.fthick = fthick;
Ti0 = (Ti-300)/F0;
Ti0 = Ti0(:,1:length(expdet.xsamp));

