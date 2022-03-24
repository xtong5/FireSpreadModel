classdef geom < handle
% geometry and visualization routines


properties
L        % length of a pixel in meters
N        % number of grid points
cx       % dimensionless x-coordinates of grid points
cy       % dimensionless y-coordinates of grid points
dx       % dimensionless width of a single cell

% state = 0 => Cell has no fuel (either burnt or never existed)
% state = 1 => Cell is not on fire but has reamining fuel
% state = 2 => Cell is currently on fire
state    % state of the burning cell
burnTime % amount of time cell has been burning
flameOut % amount of time before cell burnsout

velx    % dimensionless x component of wind
vely    % dimensionless y component of wind

end % properties

methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = geom(prams) %initial properties as structure

o.N = prams.N;
o.L = prams.L;
o.state = ones(o.N);
o.burnTime = zeros(o.N);
o.flameOut = prams.flameOut;

[o.cy,o.cx] = meshgrid(1:o.N,1:o.N);
o.dx = prams.dx;

% initialize velocity to be all zero
o.velx = zeros(o.N);
o.vely = zeros(o.N);

end % geom: constructor


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initialState(o,IC)

if strcmp(IC,'head')
  HalfWidth = floor(o.N/6);
  o.state(ceil(o.N/2)-HalfWidth:ceil(o.N/2)+HalfWidth,ceil(o.N/6)) = 2;
elseif strcmp(IC,'back')
  HalfWidth = floor(o.N/4);
  o.state(ceil(o.N/2)-HalfWidth:ceil(o.N/2)+HalfWidth,ceil(5*o.N/6)) = 2;
elseif strcmp(IC,'broken')
  HalfWidth = floor(o.N/4.2);
  o.state(ceil(o.N/4) - HalfWidth:ceil(o.N/4)+HalfWidth,ceil(o.N/4)) = 2;
  o.state(ceil(3*o.N/4) - HalfWidth:ceil(3*o.N/4)+HalfWidth,ceil(o.N/4)) = 2;
elseif strcmp(IC,'spot')
%   o.state(ceil(o.N/2),ceil(o.N/2)) = 2;
  o.state(ceil(o.N*1/5),ceil(o.N*3/4)) = 2;
  o.state(ceil(o.N*2/5),ceil(o.N*3/4)) = 2;
  o.state(ceil(o.N*3/5),ceil(o.N*3/4)) = 2;
  o.state(ceil(o.N*4/5),ceil(o.N*3/4)) = 2;
  o.state(ceil(o.N*1/4),ceil(o.N/2)) = 2;
  o.state(ceil(o.N*2/4),ceil(o.N/2)) = 2;
  o.state(ceil(o.N*3/4),ceil(o.N/2)) = 2;
  o.state(ceil(o.N*1/3),ceil(o.N*1/4)) = 2;
  o.state(ceil(o.N*2/3),ceil(o.N*1/4)) = 2;
elseif strcmp(IC,'drip')
  nignite = ceil(o.N/2);
  o.state(randi([1 o.N],nignite,1),ceil(o.N/6)) = 2;
elseif strcmp(IC,'flanks')
  HalfWidth = floor(o.N/4.2);
  o.state(ceil(o.N/2)-(4:6),...
      ceil(o.N/4) - HalfWidth:ceil(o.N/4)+HalfWidth) = 2;
  o.state(ceil(o.N/2)+(4:6),...
      ceil(o.N/4) - HalfWidth:ceil(o.N/4)+HalfWidth) = 2;
elseif strcmp(IC,'striphead')
  HalfWidth = floor(o.N/4);
  o.state(ceil(o.N/2)-HalfWidth:ceil(o.N/2)+HalfWidth,ceil(o.N/2)) = 2;
  o.state(ceil(o.N/2)-HalfWidth:ceil(o.N/2)+HalfWidth,ceil(o.N/2)-50) = 2;
elseif strcmp(IC,'ring')
  xc = ceil(o.N/2);
  yc = ceil(o.N/2);
  [xx,yy] = meshgrid(1:o.N,1:o.N);
  dist2 = (xx - xc).^2 + (yy - yc).^2;
  ind = find(dist2 > (0.98*o.N/6)^2 & dist2 < (1/0.98*o.N/6)^2);
  ind = datasample(ind,round(0.4*numel(ind)),'Replace',false);
  o.state(ind) = 2;
elseif strcmp(IC,'junction1')
  ind1 = sub2ind(size(o.state),...
    ceil(o.N/4):ceil(o.N/2),ceil(o.N/4):ceil(o.N/2));
  ind2 = sub2ind(size(o.state),...
    ceil(o.N/2):ceil(3*o.N/4),ceil(o.N/2):-1:ceil(o.N/4));
  o.state(ind1) = 2;
  o.state(ind2) = 2;
elseif strcmp(IC,'junction2')
  ind1 = sub2ind(size(o.state),...
    [ceil(o.N/4):ceil(o.N/2) ceil(o.N/4):ceil(o.N/2)],...
    [ceil(o.N/2):-1:ceil(o.N/4) ceil(o.N/2)-1:-1:ceil(o.N/4)-1]);
  ind2 = sub2ind(size(o.state),...
    [ceil(o.N/2):ceil(3*o.N/4) ceil(o.N/2):ceil(3*o.N/4)],...
    [ceil(o.N/4):ceil(o.N/2) ceil(o.N/4)-1:ceil(o.N/2)-1]);
  o.state(ind1) = 2;
  o.state(ind2) = 2;
elseif strcmp(IC,'junction3')
  [indx1,indy1] = o.bham(o.cx(ceil(0.2*o.N),ceil(0.2*o.N)), ...
                         o.cy(ceil(0.2*o.N),ceil(0.2*o.N)), ...
        ceil(0.6*o.N),ceil(o.N/8));
  indx1 = [indx1;indx1];
  indy1 = [indy1;indy1-1];
  ind1 = sub2ind(size(o.state),indx1,indy1);
  [indx2,indy2] = o.bham(o.cx(indx1(1),indy1(1)), ...
                         o.cy(indx1(1),indy1(1)), ...
        ceil(0.3*o.N),ceil(0.2*o.N));
  indx2 = [indx2;indx2];
  indy2 = [indy2;indy2+1];
  ind2 = sub2ind(size(o.state),indx2,indy2);
  
  o.state(ind1) = 2;
  o.state(ind2) = 2;
elseif strcmp(IC,'junction30')
  [indx1,indy1] = o.bham(100,40,10,37);
  indx1 = [indx1;indx1];
  indy1 = [indy1;indy1-1];
  ind1 = sub2ind(size(o.state),indx1,indy1);
  [indx2,indy2] = o.bham(100,40,-10,37);
  indx2 = [indx2;indx2];
  indy2 = [indy2;indy2+1];
  ind2 = sub2ind(size(o.state),indx2,indy2);
  
  o.state(ind1) = 2;
  o.state(ind2) = 2;
elseif strcmp(IC,'junction60')
  [indx1,indy1] = o.bham(100,40,21,37);
  indx1 = [indx1;indx1];
  indy1 = [indy1;indy1-1];
  ind1 = sub2ind(size(o.state),indx1,indy1);
  [indx2,indy2] = o.bham(100,40,-21,37);
  indx2 = [indx2;indx2];
  indy2 = [indy2;indy2+1];
  ind2 = sub2ind(size(o.state),indx2,indy2);
  
  o.state(ind1) = 2;
  o.state(ind2) = 2;
elseif strcmp(IC,'junction90')
  [indx1,indy1] = o.bham(100,40,37,37);
  indx1 = [indx1;indx1];
  indy1 = [indy1;indy1-1];
  ind1 = sub2ind(size(o.state),indx1,indy1);
  [indx2,indy2] = o.bham(100,40,-37,37);
  indx2 = [indx2;indx2];
  indy2 = [indy2;indy2+1];
  ind2 = sub2ind(size(o.state),indx2,indy2);
  
  o.state(ind1) = 2;
  o.state(ind2) = 2;
elseif strcmp(IC,'perturbedflank')
  ind = sub2ind(size(o.state),...
    ceil(o.N/2)+ceil(3*sin(0.3*(ceil(o.N/4):ceil(o.N/2)))),...
    ceil(o.N/4):ceil(o.N/2));
  o.state(ind) = 2;
end

end % initialState

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function blackout(o) %boundary?

o.state(1:ceil(o.N/10),:) = 0;        % left side
o.state(:,1:ceil(o.N/10)) = 0;        % bottom side
o.state(o.N-ceil(o.N/10):o.N,:) = 0;  % right side
o.state(:,o.N-ceil(o.N/10):o.N) = 0;  % top side


end % blackout


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,y] = bham(o,x0,y0,dx,dy)
% [x,y] = bham(x0,y0,dx,dy) computes the (x,y) coordinates of points
% starting at (x0,y0) and moving to the displacement vector (dx,dy)
% This function uses geometric reflections to address the different 8
% octants that can be in the direction (dx,dy)

% If slope is less than 1
if abs(dy) < abs(dx)
  % if going in positive x direction
  if dx > 0 
    [x,y] = o.bhamLine(dx,dy);
  % negate dx so that it falls in the previous case
  else
    [x,y] = o.bhamLine(-dx,dy);
    x = -x;
  end
% If slope is greater than 1
else
  if dy > 0 
    [y,x] = o.bhamLine(dy,dx);
  else
    [y,x] = o.bhamLine(-dy,dx);
    y = -y;
  end
end

% add on initial condition
x = x0 + x;
y = y0 + y;

end % bham

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xs,ys] = bhamLine(o,dx,dy)

% want to save all positions
xs = [];
ys = [];

x0 = 0; y0 = 0;
x1 = x0 + dx;
y1 = y0 + dy;

yi = 1;
if dy < 0
  yi = -1;
  dy = -dy;
end
D = 2*dy - dx;
y = y0;

xi = 1;
if dx < 0
  xi = -1;
end
for x = x0:xi:x1
  xs = [xs;x];
  ys = [ys;y];
  if D > 0
    y = y + yi;
    D = D + 2 * (dy - dx);
  else
    D = D + 2*dy;
  end
end


end % bhamLine



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vis(cells,time,k,streams,xstart,ystart) %plot

subplots = false;

if nargin == 2
  k = 1;
  streams = false;
end

indicator = zeros(size(cells.cx));
indicator(cells.state == 2) = 1;
indicator(cells.state == 1) = 2;

figure(k); clf; hold on
if subplots
  subplot(1,2,1)
end
titleStr = ['Time = ' num2str(time,'%4.2e')];
title(titleStr,'fontsize',20)
surface(cells.dx*(cells.cx-0.5)*cells.L,cells.dx*(cells.cy-0.5)*cells.L,indicator);
view(2);
shading flat;
colormap([[0 0 0];[1 0 0];[0.6 1 0.6]])
axis equal
axis([0 cells.dx*cells.N*cells.L 0 cells.dx*cells.N*cells.L])
set(gca,'fontsize',15)
xlabel('meters','fontsize',16)
ylabel('meters','fontsize',16)
xticks([0 50 100 150 200])
yticks([0 50 100 150 200])
if subplots
  subplot(1,2,2)
  h = quiver(cells.dx*(cells.cx-0.5)*cells.L,...
             cells.dx*(cells.cy-0.5)*cells.L,...
             cells.velx,cells.vely);
  set(h,'Linewidth',1)
  set(h,'color','black')
  axis equal
  axis([0 cells.dx*cells.N*cells.L 0 cells.dx*cells.N*cells.L])
  set(gca,'fontsize',15)
  xlabel('meters','fontsize',16)
  ylabel('meters','fontsize',16)
end

if streams
  if subplots
  subplot(1,2,1);
  end
  h = streamline(cells.cy,cells.cx,cells.vely,cells.velx,...
        ystart,xstart);
  for k=1:numel(h)
    XData = h(k).XData;
    YData = h(k).YData;
    h(k).XData = YData;
    h(k).YData = XData;
    h(k).ZData = 10*ones(size(XData));
    h(k).Color = [0 0 1];
  end
  set(h,'LineWidth',0.2)

  if subplots
  subplot(1,2,2);
  h = streamline(cells.cy,cells.cx,cells.vely,cells.velx,...
        ystart,xstart);
  for k=1:numel(h)
    XData = h(k).XData;
    YData = h(k).YData;
    h(k).XData = YData;
    h(k).YData = XData;
    h(k).Color = [0 0 1];
  end
  end

end

axis equal
axis([0 cells.dx*cells.N*cells.L 0 cells.dx*cells.N*cells.L])
set(gca,'fontsize',15)
xlabel('meters','fontsize',16)
ylabel('meters','fontsize',16)

end % vis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fat = fat(cells,state,dt) %first arrival time/state change?

fat = NaN*ones(cells.N);

[n,m] = find(state(:,:,1) == 2);
for k = 1:numel(n)
  fat(n(k),m(k)) = 0;
end

for ell = 2:size(state,3)
  dstate = state(:,:,ell) - state(:,:,ell-1); %from 
  [n,m] = find(dstate == 1);
  for k = 1:numel(n)
    fat(n(k),m(k)) = ell*dt;
  end
end

%fat = conv2(fat,ones(3)/9,'same');

% figure(2); clf;
% titleStr = ['First Arrival Time'];
% title(titleStr,'fontsize',20)
% surface(cells.dx*(cells.cx-0.5)*cells.L,...
%         cells.dx*(cells.cy-0.5)*cells.L,(fat));
% view(2);
% shading flat;
% colorbar
% cmap = buildcmap('kbryw');
% colormap(cmap)
% axis equal
% axis([0 cells.dx*cells.N*cells.L 0 cells.dx*cells.N*cells.L])
% set(gca,'fontsize',15)
% xlabel('meters','fontsize',16)
% ylabel('meters','fontsize',16)


end % fat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,y,fatVec] = fatVec(cells,fat,x0,y0,x1,y1)

[indx,indy] = cells.bham(x0,y0,x1-x0,y1-y0);

x = zeros(size(indx));
y = zeros(size(indx));
fatVec = zeros(size(indx));
for k = 1:numel(fatVec)
  x(k) = cells.cx(indx(k),indy(k)); 
  y(k) = cells.cy(indx(k),indy(k));
  fatVec(k) = fat(indx(k),indy(k));
end

end % fatVec

end % methods

end % classdef
