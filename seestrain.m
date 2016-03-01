function seestrain(expdet, strainstuff)

fs = 20; % fontsize
lw = 2; % linewidth

%ymax = 30e-3;
%ymax = 4e-3;
%figure
%hold all

xcap = expdet.xcap;
xsamp = expdet.xsamp;
xsub = expdet.xsub;
x = expdet.x;
strain = strainstuff.strain;
straine = strainstuff.straine;
t = strainstuff.t;

ymax = max(max(strain));

capind = 1:length(xcap);
sampind = length(xcap)+1:length(xcap) + length(xsamp);
subind = sampind(end)+1:length(x);

dxsamp = xsamp(3) - xsamp(2);
dxsub = xsub(3) - xsub(2);

% 3*beta*B/(rho*v2)

record_flag = 0; % set to 0 for NOT recording a video and 1 for recording

if record_flag == 1
    writerObj = VideoWriter('strain_simulation_filmFe_12ps_nocap.avi');
    writerObj.FrameRate = 10;
    writerObj.Quality = 100;
    open(writerObj);
    
end

%figure();
xlabel('film depth (nm)')
ylabel('change in lattice spacing')

i = 294;

for i = 1:5:length(t)-2
    %hold all
    plot(xcap, strain(i,capind))
    hold on
    plot(xsamp, strain(i,sampind), 'b', 'linewidth', lw)
    plot(xsub, strain(i,subind), 'linewidth', lw)
    hold off
    set(gca, 'fontsize', fs)
    text(0.8*x(end), 0.8*ymax, [num2str(t(i)) ' ps']);
    xlabel('film depth (nm)')
    ylabel('fractional change in lattice spacing')
    %hold off
    ylim([-ymax ymax]);
    drawnow
    pause(0.01);
    
    if record_flag == 1
        frame = getframe;
        writeVideo(writerObj,frame);
    end
    
end


if record_flag == 1
    close(writerObj);
end


end










