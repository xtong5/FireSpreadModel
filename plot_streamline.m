clear
close all

filename = 'ember_wash_withVort_MoreStream1_Sec300';
Setname = sprintf('DataSets/dataset_%s.mat',filename);
load(Setname)

save_plot_inter_streams = 1;
save_plot_final_streams = 1;

%% stream lines start location
xstart1 = linspace(22,179,50);
ystart1 = 22*ones(size(xstart1));

ystart2 = linspace(22,179,20);
xstart2 = 22*ones(size(ystart2));

ystart3 = linspace(22,179,20);
xstart3 = 179*ones(size(ystart3));

xstart = [xstart1 xstart2 xstart3];
ystart = [ystart1 ystart2 ystart3];

%% streamlines at a specified time/step
i = 10; %time step
time = prams.dt*(i-1); 
fprintf('Stream lines at %g step, %f seconds.\n', i, time)

indicator = zeros(size(cx));
indicator(state(:,:,i) == 2) = 1;
indicator(state(:,:,i) == 1) = 2;

figure(1); clf; hold on
titleStr = ['Time = ' num2str(time,'%4.2e')];
title(titleStr,'fontsize',20)
surface(prams.dx*(cx-0.5)*prams.L,prams.dx*(cy-0.5)*prams.L,indicator);
view(2);
shading flat;
colormap([[0 0 0];[1 0 0];[0.6 1 0.6]])
axis equal
axis([0 prams.dx*prams.N*prams.L 0 prams.dx*prams.N*prams.L])
set(gca,'fontsize',15)
xlabel('meters','fontsize',16)
ylabel('meters','fontsize',16)
xticks([0 50 100 150 200])
yticks([0 50 100 150 200])

h = streamline(cy,cx,vely(:,:,i),velx(:,:,i),ystart,xstart);
for k=1:numel(h)
    XData = h(k).XData;
    YData = h(k).YData;
    h(k).XData = YData;
    h(k).YData = XData;
    h(k).ZData = 10*ones(size(XData));
    h(k).Color = [0 0 1];
end
set(h,'LineWidth',0.2)

if save_plot_inter_streams
  name = sprintf('%s_streamlines_at_step%g',filename,i);
  saveas(1,name,'png')
end

% streamlines at final step/state
i = size(state,3);
time = prams.dt*(i-1); 
fprintf('Stream lines at %g (final)step, %f seconds.\n', i, time)

indicator = zeros(size(cx));
indicator(state(:,:,i) == 2) = 1;
indicator(state(:,:,i) == 1) = 2;

figure(2); clf; hold on
titleStr = ['Time = ' num2str(time,'%4.2e')];
title(titleStr,'fontsize',20)
surface(prams.dx*(cx-0.5)*prams.L,prams.dx*(cy-0.5)*prams.L,indicator);
view(2);
shading flat;
colormap([[0 0 0];[1 0 0];[0.6 1 0.6]])
axis equal
axis([0 prams.dx*prams.N*prams.L 0 prams.dx*prams.N*prams.L])
set(gca,'fontsize',15)
xlabel('meters','fontsize',16)
ylabel('meters','fontsize',16)
xticks([0 50 100 150 200])
yticks([0 50 100 150 200])

h = streamline(cy,cx,vely(:,:,i),velx(:,:,i),ystart,xstart);
for k=1:numel(h)
    XData = h(k).XData;
    YData = h(k).YData;
    h(k).XData = YData;
    h(k).YData = XData;
    h(k).ZData = 10*ones(size(XData));
    h(k).Color = [0 0 1];
end
set(h,'LineWidth',0.2)

if save_plot_final_streams
  name = sprintf('%s_streamlines_final_step%g',filename,i);
  saveas(2,name,'png')
end