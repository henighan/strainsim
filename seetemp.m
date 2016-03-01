function seetemp(expdet, Te, Ti)

lw = 2; % linewidths for plots
fs = 20; % fontsize for plots


dxsamp = expdet.dxsamp;
xsamp = dxsamp:dxsamp:dxsamp*size(Te, 2);

t = expdet.t;

ymax = max(Te(:));
ymin = min(Te(:));

figure(1)
title('blue is lattice, red is electrons')
for i = 1:5:length(t)-2
    %hold all
    plot(xsamp, Te(i,:), 'r', 'linewidth', lw)
    hold on
    plot(xsamp, Ti(i,:), 'b', 'linewidth', lw)
    %hold off
    text(0.8*xsamp(end), 0.8*ymax, [num2str(t(i)) ' ps'], ...
        'fontsize', fs);
    xlabel('film depth (nm)', 'fontsize', fs)
    ylabel('Temperature (K)', 'fontsize', fs)
    set(gca, 'fontsize', fs)
    hold off
    ylim([0.8*ymin ymax]);
    drawnow
    pause(0.01);
    
    
end


%%
%{
figure
hold all
lw = 2;
fs = 20;
plot(expdet.t, mean(Te, 2), 'linewidth', lw)
plot(expdet.t, mean(Ti, 2), 'linewidth', lw)
title('Depth-averaged temperatures')
%plot(expdet.t, Te(:,1), 'linewidth', lw)
%plot(expdet.t, Ti(:,1), 'linewidth', lw)
%title('Surface temperatures')
%plot(t, Ts(:,1))
xlabel('time (ps)')
ylabel('temperature (Kelvin)')
legend('Electrons', 'Lattice')
set(gca, 'fontsize', fs)
%}
%%

end
