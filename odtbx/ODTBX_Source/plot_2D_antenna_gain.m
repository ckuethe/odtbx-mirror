function plot_2D_antenna_gain(fname)

% fname = 'S41P1L1Rpublic.txt';
gain = load(fname);

az = gain(1,2:end);
el = gain(2:end,1);
gain = gain(2:end,2:end);

% repeat the 0 deg column so figure wraps around
gain = [gain gain(:,1)];
az = [az 360];


figure;
[xi,yi,zi]=polarplot3d(gain,'PlotType','surf',...
    'AngularRange',az*pi/180,'RadialRange',el,...
    'TickSpacing',15,'PolarDirection','ccw');
cax = caxis;
caxis([-30 15]);
colormap jet
h=colorbar;
set(get(h,'ylabel'),'String','Gain (dB)');
shading INTERP
zlabel('Gain (dB)')
set(gca,'Xtick',[-90:15:90],'Ytick',[-90:15:90]);
axis square
xlabel('Off-Boresight Angle (deg)')
title(fname);
figname = char(strcat(fname,'.jpg'));
% print('-djpeg90', figname)
set(gca,'CameraPosition',[0 0 100]);  % Look along the z-axis for 2D




end