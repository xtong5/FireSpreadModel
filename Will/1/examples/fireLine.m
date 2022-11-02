clear

addpath ../src

usePlot = true;
savePlot = false;
saveFrames = false;
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

prams.uwind = 0.1; % dimensionless speed of wind from the south
prams.windDir = 45; % wind direction relative to north in degrees
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
prams.probIgnite = 4e-1;  

% create object to store information (state, burn time, flame out,
% velocity) of each cell
cells = geom(prams);

% create ignition pattern
cells.initialState('spot');
%cells.initialState('user',flameOutTime);
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
  
  % calculate vertical velocity
  N = prams.N;
  velz = zeros(N, N);
  for x=2:N-1
   for y=2:N-1
    %velz(x,y) = cells.velx(x+1,y) - cells.velx(x-1,y) + cells.vely(x,y+1) - cells.velx(x,y-1);
    velz(x,y) = cells.velx(x+1,y) - cells.velx(x-1,y) + cells.vely(x,y+1) - cells.velx(x,y-1);
   end
  end
  %velz = (-1/2) * velz;
  velz = (-1/2) * velz;

  % display graph of vertical velocity
  % correct for plotting issues
  velz2 = rot90(fliplr(velz));
  figure(2); clf; hold on
  subplot(2,2,1);
  imagesc(velz2);
  title("Vertical Velocity");
  set(gca, "YDir", "normal");
  axis equal tight;
  colorbar;

  % calculate vorticity
  vortx = zeros(N, N);
  vorty = zeros(N, N);
  vortz = zeros(N, N);
  for x=2:N-1
   for y=2:N-1
    vortx(x,y) = (velz(x,y+1) - velz(x,y-1))/2 - cells.vely(x,y);
    vorty(x,y) = cells.velx(x,y) - (velz(x+1,y) - velz(x-1,y))/2;
    vortz(x,y) = (1/2) * (cells.vely(x+1,y) - cells.vely(x-1,y) - cells.velx(x,y+1) + cells.velx(x,y-1));
   end
  end
  % not sure why, but vortz appears to have two "layers" of erroneous data,
  % instead of the single layer found in x and y
  vortz(:,2) = 0;
  vortz(:,end-1) = 0;
  vortz(2,:) = 0;
  vortz(end-1,:) = 0;

  vortxy = reshape([vortx, vorty], N, N, 2);
  horizVort = vecnorm(vortxy, 2, 3);

  % display vorticity
  % correct for plotting issues
  % vortx2 = rot90(fliplr(vortx));
  % vorty2 = rot90(fliplr(vorty));
  horizVort2 = rot90(fliplr(horizVort));
  vortz2 = rot90(fliplr(vortz));
  % subplot(2,2,2);
  % imagesc(vortx2);
  % set(gca, "YDir", "normal");
  % axis equal tight;
  % colorbar;
  % subplot(2,2,3);
  % imagesc(vorty2);
  % set(gca, "YDir", "normal");
  % axis equal tight;
  % colorbar;
  subplot(2,2,2);
  imagesc(horizVort2);
  title("Magnitude of Horizontal Vorticity");
  set(gca, "YDir", "normal");
  axis equal tight;
  colorbar;
  subplot(2,2,3);
  imagesc(vortz2);
  title("Vertical Vorticity");
  set(gca, "YDir", "normal");
  axis equal tight;
  colorbar;

  % calculate helicity
  vel = reshape([cells.velx cells.vely velz], N, N, 3);
  vort = reshape([vortx vorty vortz], N, N, 3);
  H = dot(vel, vort, 3);

  % display helicity
  % correct for plotting issues
  H2 = rot90(fliplr(H));
  subplot(2,2,4)
  imagesc(H2);
  title("Helicity")
  set(gca, "YDir", "normal");
  axis equal tight;
  colorbar;
  
  % display velocity
%   subplot(2,2,4)
%   vel2 = rot90(fliplr(vel(:,:,1:2)));
%   imagesc(vecnorm(vel2(:,:,1:2), 2, 3));
%   set(gca, "YDir", "normal");
%   axis equal tight;
%   colorbar;

  state(:,:,i) = cells.state;
  velx(:,:,i) = cells.velx;
  vely(:,:,i) = cells.vely;
  i = i + 1;

  % visualize the current state
  if usePlot
    cells.vis(time,1,streams,xstart,ystart,flameOutTime);
    if saveFrames
      name = sprintf('frames/frame_%s_%03i',filename,i-1);
      saveas(1,name,'png');
      saveas(2,strcat(name,'_fig2'),'png');
    end
  end
  % update time
  time = time + prams.dt;
  % pause(.01)

  % update the state of the cells and the amount of time they've been
  % burning
  os.updateState(cells);
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
