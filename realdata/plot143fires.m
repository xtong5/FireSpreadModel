clear

load('selectedfires.mat')

N = size(T,3);


% filename = sprintf('143fires');
% Vname = sprintf('%s',filename);

%     v = VideoWriter(Vname,'MPEG-4');
%     v.FrameRate = 2;
%     open(v)
    
for k = 1:N

% figure(1);clf;
t = T(:,1,k);
A = T(:,2,k);
titleStr = sprintf('Area vs Time, fire=%s',Names(k));
title(titleStr,'fontsize',20)
loglog(t,A,'-','LineWidth',2)
hold on

%     frame = getframe(gcf);
%     writeVideo(v,frame);
%  


end


% 
% ylim([0, 0.9])
xlim([1e-3, 1e-2])
xlabel('Time scaled','fontsize',16)
ylabel('Area scaled','fontsize',16)


%    close(v)



