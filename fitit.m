function fitit(j, nsteps)

cd('~/strainsim2/');

%nsteps = 10;
ub = 5*3.6;
stepsize = ub/(nsteps);

elclist =0:stepsize:ub; 

elstrcoeff = elclist(j);

disp(['elstrcoeff is ' num2str(elstrcoeff)])

%matlabpool('open', 4)

%load('temps15f23nm');

geometry = load_geometry('CSPAD', 'Fe');

expdet = load_expdet();

F0 = 15;

[Ti, Te, Ts] = spinstempsim2(expdet.samp, expdet, F0);


fthick = 22.6;
[strainstuff, expdet] = strainfunc_vthick(expdet, Te, Ti, fthick);
strainstuff.fthick = fthick;
Ti0 = (Ti-300)/F0;
Ti0 = Ti0(:,1:length(expdet.xsamp));

F0 = 15;
%elstrcoeff = 0;
atl = 660;
phioffset = -0.66;
ksbw0 = 0.006;

%paramvec0 = [F0 elstrcoeff phioffset atl];
paramvec0 = [F0 phioffset atl];

chi2handle = @biggerchi2;
paramvecout = paramvec0;
options = optimset('MaxFunEvals',1e10, 'MaxIter', 1e10);
%
%tic
paramvecout = ...
    fminsearch(@(paramvec) biggerchi2(strainstuff, geometry, Ti0, elstrcoeff, paramvec), ...
    paramvec0, options)
%toc

%}

[chi2, t0list] = biggerchi2(strainstuff, geometry, Ti0, elstrcoeff, paramvecout);

fitstuff = struct(  'paramvecout',paramvecout,...
                    'chi2',chi2,...
                    'elstrcoeff',elstrcoeff,...
                    't0list',t0list);
                
save(['~/strainsim2/fitoutput/fitstuff' num2str(j) '.mat'], '-struct', 'fitstuff');

end
