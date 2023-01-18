% clear

addpath ../src

usePlot = true;
savePlot = false;
saveFrames = true;
saveData = false;

filename = 'fireLine';
Vname = sprintf('test_%s',filename);

% velocity scale is 1m/s
%prams.V = 1;
% length scale is 1m (typical size of a single pixel)
prams.L = 1;

% number of points in x and y direction. Also sets the domain size
prams.N = 201; 
prams.dx = 1; % dimensionless size of a cell

prams.uwind = 0; % dimensionless speed of wind from the south
prams.windDir = 310; % wind direction relative to north in degrees
prams.s = prams.uwind/4; 
% dimensionless speed of the fire spread rate due only to wind

% 2x uwind produces lateral spread effects strength of omega in Sharples
% paper = 0.3 (scaling?)
prams.vortStrength = 0;
% prams.vortStrength = 0.5*prams.uwind/prams.dx^2;
% prams.vortStrength = 0.5*prams.uwind/prams.dx^2; % strength of the vorticity
% prams.vortStrength = -2*prams.s/prams.dx^2; % strength of the vorticity

% 0.9 base runs %4. %4. %0.3; 1x uwind has large effect. Adding 1 causes
% blow-up at uwind = 5 strength of nu in Sharples paper. The minus sign
% to make it a sink is accounted for later in the code
%prams.sinkStrength = (0.2*(exp(-prams.uwind) + 0.8*prams.uwind))/...
%  prams.dx^2;
% prams.sinkStrength = 0.0;
prams.sinkStrength = 0.2; % for uwind = 1
% prams.sinkStrength = 0.23; % for uwind = 1
% prams.sinkStrength = 0.35; % for uwind = 2

% time required for cell to travel 3 grid points
if prams.s ~= 0
  prams.dt = 3*prams.dx/prams.s; 
else
  prams.dt = 2.0;
end
prams.T = 10000; % time horizon
prams.ntime = ceil(prams.T/prams.dt);

prams.spotting = false;
prams.emberWash = false;


if 0
load('fuelmaps/fuels_30-70split_eps1em2.mat');
end
if 0
load('fuelmaps/fuels_30-70split_eps1em3.mat');
end
if 0
load('fuelmaps/fuels_30-70split_eps1em4.mat');
end
if 0
load('fuelmaps/fuels_50-50split_eps1em2.mat');
end
if 0
load('fuelmaps/fuels_50-50split_eps1em3.mat');
end
if 0
load('fuelmaps/fuels_50-50split_eps1em4.mat');
end
if 1
load('fuelmaps/fuels_uniform.mat');
end
fuels = sign(fuels);
fuels = fuels(1:prams.N,1:prams.N);
% max and minimum burn times in seconds
maxBurn = 60;
minBurn = 15;
flameOutTime = zeros(prams.N);
flameOutTime(fuels == 1) = maxBurn;
flameOutTime(fuels == -1) = minBurn;


% Time steps that cell can burn for ...in zero windXX make ~ 1/wind
% number of dimensionless time steps before a cell burns out
prams.flameOut = ceil(flameOutTime/prams.dt); 

% probability that a cell ignites a neighbor to model diffusion.
% Probabilities of igniting diagonal neighbors is scaled by sqrt(2) to
% account for additional distance. Ignition due to diffusion only occurs
% just before a cell burns out
% prams.probIgnite = 4e-1;
prams.probIgnite = 0.4;

% set logistic growth rate of uphill diffusion weight
% ignition probability is weighted to be more likely uphill, and scales
% on a logistic curve with slope. 
if exist('strength', 'var') % allows running from command line or caller.m
  prams.strengthSlopeDiffuse = strength;
else
  prams.strengthSlopeDiffuse = 50;
end

% create object to store information (state, burn time, flame out,
% velocity) of each cell
cells = geom(prams);

% create ignition pattern
cells.initialState('head');
% cells.initialState('user',flameOutTime);
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
xstartBot = linspace(6,prams.N-5,50);
ystartBot = 6*ones(size(xstartBot));

ystartLeft = linspace(6,prams.N-5,50);
xstartLeft = 6*ones(size(ystartLeft));

ystartRight = linspace(6,prams.N-5,50);
xstartRight = (prams.N-5)*ones(size(ystartRight));

xstart = [xstartBot xstartLeft xstartRight];
ystart = [ystartBot ystartLeft ystartRight];
%xstart = [xstartBot xstartLeft];
%ystart = [ystartBot ystartLeft];

% select the topographic map to use
if 0
  % generate pseudorandom terrain
  topo = 1-randTerrain(5, prams.N);
end
if 0
  % linear ridge topography, vertical
  load('linear_ridge.mat');
end
if 1
  % linear ridge topography, horizontal
  load('linear_ridge.mat');  
  topo = topo';
end
topoPlot = rot90(fliplr(topo)); % correct for plotting issues

% calculating the gradient from topography
grady = zeros(prams.N);
gradx = zeros(prams.N);
grady(2:end-1,2:end-1) = topo(2:end-1,3:end) - topo(2:end-1,1:end-2);
gradx(2:end-1,2:end-1) = topo(3:end,2:end-1) - topo(1:end-2,2:end-1);
grady = grady .* 0.5;
gradx = gradx .* 0.5;
grad = reshape([gradx, grady], prams.N, prams.N, 2);

i = 1;
% main time stepping loop
while time < prams.T && nnz(cells.state == 2)
  fprintf('t = %4.2e of T = %4.2e\n', time, prams.T)
  % sink and vorticity terms due to cells that are on fire. Both are
  % divided by dx to account for the fact that the (x,y) divergence is
  % dw/dz
  sinkForce = os.sinkTerm(cells);
  vortForce = os.vortTerm(cells);

  % compute velocity due to sinks
  [psi,psix,psiy] = os.PoissonSolverNeumann(sinkForce); 

  % compute velocity due to vorticity
  [psi,etax,etay] = os.PoissonSolverNeumann(+vortForce);

  % add the different terms in the velocity
  [cells.velx,cells.vely] = os.computeVelocity(...
        psix,psiy,etax,etay);

  state(:,:,i) = cells.state;
  velx(:,:,i) = cells.velx;
  vely(:,:,i) = cells.vely;
  i = i + 1;

  % visualize the current state
  if usePlot
    cells.vis(time,1,streams,xstart,ystart,flameOutTime,topoPlot,grad);
    if saveFrames
      name = sprintf('frames/frame_%s_%04i',filename,i-1);
      saveas(1,name,'png')
    end
  end
  % update time
  time = time + prams.dt;
  % pause(.01)

  % update the state of the cells and the amount of time they've been
  % burning
  os.updateState(cells, topo);
end

return

% save final state and velocities
state(:,:,prams.ntime+1) = cells.state;
velx(:,:,prams.ntime+1) = cells.velx;
vely(:,:,prams.ntime+1) = cells.vely;

% visualize the history and save the final frame if savePlot is true
if usePlot || savePlot
  cells = geom(prams);
  for k = 1:prams.ntime+1
    cells.state = state(:,:,k);
    cells.velx = velx(:,:,k);
    cells.vely = vely(:,:,k);
    time = (k-1)*prams.dt;
    cells.vis(time,1,streams,xstart,ystart,flameOutTime);
    pause(0.01)
  end
  if savePlot
    name = sprintf('%s_final_state',filename);
    saveas(1,name,'png')
  end
  %save('fireline.mat','prams','state','velx','vely');
end

if usePlot || savePlot
  fat = cells.fat(state,prams.dt);
  if savePlot
    name = sprintf('%s_FAT',filename);
    saveas(2,name,'png')
  end
end

% compute first arrival time along a particular vector
if 0
x0 = 100; y0 = 40;
x1 = 100; y1 = 80;
[x,y,fatVec] = cells.fatVec(fat,x0,y0,x1,y1);
end


%%save data
if saveData
  N_state = size(state,3);
  burnMap = NaN*ones(prams.N,prams.N,N_state);
  fuelMap = NaN*ones(prams.N,prams.N,N_state);

  for i = 1:N_state
    burning = NaN*ones(prams.N,prams.N);
    fuel = zeros(prams.N,prams.N);
    burning(state(:,:,i) == 2) = 1;
    fuel(state(:,:,i) == 2) = 1; % burning=with fuel
    fuel(state(:,:,i) == 1) = 1;
    burnMap(:,:,i) = burning;
    fuelMap(:,:,i) = fuel;
  end

  Name_data = sprintf('dataset_%s_Sec%g.mat',filename,prams.T);
  cx = cells.cx; cy = cells.cy;
  save(Name_data,'cx','cy','fat','prams','xstart','ystart','state',...
  'velx','vely','fuelMap','burnMap');
end
