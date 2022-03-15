clear

addpath src_ember_wash

usePlot = true;
saveVideo = true;

if saveVideo
    v = VideoWriter('test_ember_wash_withVort_NoEmber','MPEG-4');
    v.FrameRate = 2;
    open(v)
end


% velocity scale is 1m/s
%prams.V = 1;
% length scale is 1m (typical size of a single pixel)
prams.L = 1;

% number of points in x and y direction. Also sets the domain size
prams.N = 201; 
prams.dx = 1; % dimensionless size of a cell

prams.uwind = 4; % dimensionless speed of wind from the south
prams.s = prams.uwind/4; 
% dimensionless speed of the fire spread rate due only to wind

% 2x uwind produces lateral spread effects strength of omega in Sharples
% paper = 0.3 (scaling?) 
% prams.vortStrength = 0;
prams.vortStrength = 0.5*prams.uwind/prams.dx^2; % strength of the vorticity
%prams.vortStrength = -2*prams.s/prams.dx^2; % strength of the vorticity

% 0.9 base runs %4. %4. %0.3; 1x uwind has large effect. Adding 1 causes
% blow-up at uwind = 5 strength of nu in Sharples paper. The minus sign
% to make it a sink is accounted for later in the code
%prams.sinkStrength = (0.2*(exp(-prams.uwind) + 0.8*prams.uwind))/...
%  prams.dx^2;
prams.sinkStrength = 0.7; % for uwind = 1
%prams.sinkStrength = 0.23; % for uwind = 1
%prams.sinkStrength = 0.35; % for uwind = 2

%prams.vortStrength = 0;
%prams.sinkStrength = 0;

% time required for cell to travel 3 grid points
if prams.s ~= 0
  prams.dt = 3*prams.dx/prams.s; 
else
  prams.dt = 2.0;
%  prams.dt = 4;
end
prams.T = 300;
prams.ntime = ceil(prams.T/prams.dt);

% Time steps that cell can burn for ...in zero windXX make ~ 1/wind
% number of dimensionless time steps before a cell burns out
flameOutTime = 15; % time to flameout in seconds
prams.flameOut = ceil(flameOutTime/prams.dt); 

% probability that a cell ignites a neighbor to model diffusion.
% Probabilities of igniting diagonal neighbors is scaled by sqrt(2) to
% account for additional distance. Ignition due to diffusion only occurs
% just before a cell burns out
prams.probIgnite = 4e-1;  

% create object to store information (state, burn time, flame out,
% velocity) of each cell
cells = geom(prams);

% create ignition pattern
% cells.initialState('junction30');
cells.initialState('spot');
% state = 0 => Cell has no fuel (either burnt or never existed)
% state = 1 => Cell is not on fire but has reamining fuel
% state = 2 => Cell is currently on fire

% black out the boundary
cells.blackout;

% space to save the state at each time step
state = zeros(prams.N,prams.N,prams.ntime+1);
velx = zeros(prams.N,prams.N,prams.ntime+1);
vely = zeros(prams.N,prams.N,prams.ntime+1);

% object for the PDE solvers and the Bresenham algorithm
os = solvers(prams);

% current time
time = 0;

streams = true;
xstart = linspace(22,179,100);
ystart = 22*ones(size(xstart));

%% first step
i = 1;
disp(time)
sinkForce = os.sinkTerm(cells);
vortForce = os.vortTerm(cells);
[psi,psix,psiy] = os.PoissonSolverNeumann(sinkForce); 
[psi,etax,etay] = os.PoissonSolverNeumann(+vortForce);
cos = zeros(prams.N,prams.N); sin = cos;
[cells.velx,cells.vely] = os.computeVelocity(...
        psix,psiy,etax,etay,cos,sin);
state(:,:,i) = cells.state;
velx(:,:,i) = cells.velx;
vely(:,:,i) = cells.vely;
i = i + 1;
time = time + prams.dt;
os.updateState(cells);


while time < prams.T
  disp(time)
% If there is a prescribed number of time steps
%for i = 1:prams.ntime
  % sink and vorticity terms due to cells that are on fire. Both are
  % divided by dx to account for the fact that the (x,y) divergence is
  % dw/dz
  sinkForce = os.sinkTerm(cells);
  vortForce = os.vortTerm(cells);

  % compute velocity due to sinks
  [psi,psix,psiy] = os.PoissonSolverNeumann(sinkForce); 
%  [psi,psix,psiy] = os.PoissonSolverDirichlet(sinkForce); 

  % compute velocity due to vorticity
  [psi,etax,etay] = os.PoissonSolverNeumann(+vortForce);
%  [psi,etax,etay] = os.PoissonSolverDirichlet(+vortForce);

  % add the different terms in the velocity
  hypo = (velx(:,:,i-1).^2+vely(:,:,i-1).^2).^(0.5);
  cos = velx(:,:,i-1)./hypo; sin = vely(:,:,i-1)./hypo; 
  [cells.velx,cells.vely] = os.computeVelocity(...
        psix,psiy,etax,etay,cos,sin);

  state(:,:,i) = cells.state;
  velx(:,:,i) = cells.velx;
  vely(:,:,i) = cells.vely;
  i = i + 1;

  % visualize the current state
  if usePlot
    cells.vis(time,1,streams,xstart,ystart);
    if saveVideo
        frame = getframe(gcf);
        writeVideo(v,frame);
    end
  end
%  fprintf('time step %i! PAUSED\n',i)
  % update time
  time = time + prams.dt;
  pause(.01)
%  pause

%  figure(2);
%  os.FormStreams(cells,xstart,ystart);
%  figure(1)
%  pause

  % update the state of the cells and the amount of time they've been
  % burning
  os.updateState(cells);
end

cx = cells.cx; cy = cells.cy;
state(:,:,prams.ntime+1) = cells.state;
velx(:,:,prams.ntime+1) = cells.velx;
vely(:,:,prams.ntime+1) = cells.vely;

% visualize the final state
% if usePlot
%   cells.vis(time,1,streams,xstart,ystart);
% end
% pause(.01)

cells = geom(prams);
for k = 1:prams.ntime+1
 cells.state = state(:,:,k);
 cells.velx = velx(:,:,k);
 cells.vely = vely(:,:,k);
 time = (k-1)*prams.dt;
 cells.vis(time);
 pause(0.01)
end
%save('fireline.mat','prams','state','velx','vely');

fat = cells.fat(state,prams.dt);

% x0 = 45; y0 = 42;
% %x1 = 102; y1 = 69;
% x1 = 160; y1 = 96;

x0 = 100; y0 = 40;
x1 = 100; y1 = 80;

[x,y,fatVec] = cells.fatVec(fat,x0,y0,x1,y1);

% figure(3); clf;
% plot(y,fatVec);

% save('fireData.mat','cx','cy','fat','fatVec','prams','x','y','xstart','ystart');

figure(2); clf;
titleStr = ['First Arrival Time'];
title(titleStr,'fontsize',20)
surface(cells.dx*(cells.cx-0.5)*cells.L,...
        cells.dx*(cells.cy-0.5)*cells.L,(fat));
view(2);
shading flat;
colorbar
cmap = buildcmap('kbryw');
colormap(cmap)
axis equal
axis([0 cells.dx*cells.N*cells.L 0 cells.dx*cells.N*cells.L])
set(gca,'fontsize',15)
xlabel('meters','fontsize',16)
ylabel('meters','fontsize',16)
yticks([0 50 100 150 200])

if saveVideo
    frame = getframe(gcf);
    writeVideo(v,frame);
    close(v)
end