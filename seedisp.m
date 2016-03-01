function seedisp(strainstuff)

%ymax = 3;

figure();
xlabel('film depth (nm)')
ylabel('dispacement from equilibrium position (Angstroms)')

dxsamp = strainstuff.xsamp(2)-strainstuff.xsamp(1);
sampind = strainstuff.sampind;
strain = strainstuff.strain-strainstuff.staticstrain;
xsamp = strainstuff.xsamp;
t = strainstuff.t;

u = -10*dxsamp*fliplr(cumsum(fliplr(strain(:,sampind)), 2));
%u = -10*dxsub*fliplr(cumsum(fliplr(strain(:,subind)), 2));

ymax = max(abs(u(:)));

for i = 1:length(t)-2

    plot(xsamp, u(i,:), 'r')
    %plot(xsub, u(i,:), 'r')

    text(0.8*xsamp(end), -0.8*ymax, [num2str(t(i)) ' ps']);
    xlabel('film depth (nm)')
    ylabel('dispacement from equilibrium position (Angstroms)')
    %hold off
    ylim([-ymax ymax]);
    drawnow
    pause(0.01);
    
end

end