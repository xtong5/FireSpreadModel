classdef solvers
% different PDE solvers


properties
N
dx
dt
uwind
windDir
s
sinkStrength
vortStrength
probIgnite
L
nullVec
spotting
emberWash
strucLeft
strucBot
strucSize

end % properties

methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = solvers(prams)

o.N = prams.N;
o.dx = prams.dx;
o.dt = prams.dt;
o.uwind = prams.uwind;
o.windDir = prams.windDir*pi/180; % convert wind direction to radians
o.s = prams.s;
o.sinkStrength = prams.sinkStrength;
o.vortStrength = prams.vortStrength;
o.probIgnite = prams.probIgnite;
o.spotting = prams.spotting;
o.emberWash = prams.emberWash;

if isfield(prams,'strucLeft')
  o.strucLeft = prams.strucLeft;
  o.strucBot = prams.strucBot;
  o.strucSize = prams.strucSize;
else
  o.strucLeft = [];
  o.strucBot =  [];
  o.strucSize = [];
end

[o.L,o.nullVec] = o.PoissonSolverNeumannMatrix;

end % solvers: constructor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sink = sinkTerm(o,cells)

sink = zeros(cells.N,cells.N);
sink(cells.state == 2) = -o.sinkStrength;

end % sinkTerm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vort = vortTerm(o,grid)

vort = zeros(grid.N,grid.N);

% search along each row at cells that are burning. Set the point to the
% right and left to be +1 and -1 so that it forms the desired vorticity
for k = 1:grid.N
  ind = find(grid.state(:,k) == 2);
  if numel(ind) ~= 0
    % For broken fire lines, this will put a positive vorticity just to
    % the left and a negative vorticity just to the right of each
    % connected component. There is no vorticity force inside the
    % fireline
    vort(ind-1,k) = vort(ind-1,k) + o.vortStrength;
    vort(ind+1,k) = vort(ind+1,k) - o.vortStrength;
    vort(ind,k) = 0;
  end
end

% account for the fact that dw/dz depends on the grid scale (dz)
%vort = vort/o.dx;

end % vortTerm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [L,nullVec] = PoissonSolverNeumannMatrix(o)
% Compupte the matrix for 2D Poisson equation with Neumann boundary
% conditions

N = o.N;

% II needed for top and bottom boundary conditions
II = speye(N);
% IIt needed for interaction between points above and below. It is zero
% in the first and last rows because we are enforcing no flux in the
% horizontal (normal) direction at these points
IIt = II; IIt(1,1) = 0; IIt(end,end) = 0;

% D is the matrix that handles diagonal terms and neighbors to the left
% and right. First and last rows are zero since these are the 
e = ones(N,1);
D = spdiags([e -4*e e],(-1:1),N,N);
D(1,1) = -1*o.dx; D(1,2) = 1*o.dx;
D(N,N-1) = 1*o.dx; D(N,N) = -1*o.dx;

L = spalloc(N^2,N^2,6*N);
% Neumann boundary condition at bottom
L(1:N,1:N) = -II/o.dx;
L(1:N,N+1:2*N) = II/o.dx;
% Neumann boundary condition at top
L(N*(N-1)+1:end,N*(N-1)+1:end) = -II/o.dx;
L(N*(N-1)+1:end,N*(N-2)+1:N*(N-1)) = II/o.dx;

% points below
for k = 2:N-1
  icol = (k-2)*N+1;
  irow = (k-1)*N+1;
  L(irow:irow+N-1,icol:icol+N-1) = IIt/o.dx^2;
end

% points above
for k = 2:N-1
  icol = k*N+1;
  irow = (k-1)*N+1;
  L(irow:irow+N-1,icol:icol+N-1) = IIt/o.dx^2;
end

% diagonal and points to left and right
for k = 2:N-1
  icol = (k-1)*N+1;
  irow = (k-1)*N+1;
  L(irow:irow+N-1,icol:icol+N-1) = D/o.dx^2;
end

[cy,cx] = meshgrid(1:o.N,1:o.N);

eta = 1; % penalty parameter
% modifications to L due to the obstacles
indLeft = o.strucLeft;
indBot = o.strucBot;
len = o.strucSize;
dx = o.dx;
for j = 1:numel(indLeft)
  % normal derivative at lower left corner
  ind1 = sub2ind([N,N],indLeft(j),indBot(j));
  ind2 = sub2ind([N,N],indLeft(j),indBot(j)-1);
  L(ind1,:) = 0;
  L(ind1,ind1) = 0*L(ind1,ind1) - eta*1/dx;
  L(ind1,ind2) = 0*L(ind1,ind2) + eta*1/dx;

  % normal derivative at upper left corner
  ind1 = sub2ind([N,N],indLeft(j),indBot(j)+len);
  ind2 = sub2ind([N,N],indLeft(j),indBot(j)+len+1);
  L(ind1,:) = 0;
  L(ind1,ind1) = 0*L(ind1,ind1) - eta*1/dx;
  L(ind1,ind2) = 0*L(ind1,ind2) + eta*1/dx;

  % normal derivative at lower right corner
  ind1 = sub2ind([N,N],indLeft(j)+len,indBot(j));
  ind2 = sub2ind([N,N],indLeft(j)+len,indBot(j)-1);
  L(ind1,:) = 0;
  L(ind1,ind1) = 0*L(ind1,ind1) - eta*1/dx;
  L(ind1,ind2) = 0*L(ind1,ind2) + eta*1/dx;

  % normal derivative at upper right corner
  ind1 = sub2ind([N,N],indLeft(j)+len,indBot(j)+len);
  ind2 = sub2ind([N,N],indLeft(j)+len,indBot(j)+len+1);
  L(ind1,:) = 0;
  L(ind1,ind1) = 0*L(ind1,ind1) - eta*1/dx;
  L(ind1,ind2) = 0*L(ind1,ind2) + eta*1/dx;

  % start of imposing Neumann boundary condition along the sides
  % (non-corners)
  for k = indLeft(j)+1:indLeft(j) + len - 1
    % update indices for normal derivative along bottom side
    ind1 = sub2ind([N,N],k,indBot(j));
    ind2 = sub2ind([N,N],k,indBot(j)-1);
    L(ind1,:) = 0;
    L(ind1,ind1) = 0*L(ind1,ind1) - eta*1/dx;
    L(ind1,ind2) = 0*L(ind1,ind2) + eta*1/dx;

    % update indices for normal derivative along top side
    ind1 = sub2ind([N,N],k,indBot(j)+len);
    ind2 = sub2ind([N,N],k,indBot(j)+len+1);
    L(ind1,:) = 0;
    L(ind1,ind1) = 0*L(ind1,ind1) - eta*1/dx;
    L(ind1,ind2) = 0*L(ind1,ind2) + eta*1/dx;
  end

  for k = indBot(j)+1:indBot(j) + len - 1
    % update indices for normal derivative along left side
    ind1 = sub2ind([N,N],indLeft(j),k);
    ind2 = sub2ind([N,N],indLeft(j)-1,k);
    L(ind1,:) = 0;
    L(ind1,ind1) = 0*L(ind1,ind1) - eta*1/dx;
    L(ind1,ind2) = 0*L(ind1,ind2) + eta*1/dx;

    % update indices for normal derivative along right side
    ind1 = sub2ind([N,N],indLeft(j)+len,k);
    ind2 = sub2ind([N,N],indLeft(j)+len+1,k);
    L(ind1,:) = 0;
    L(ind1,ind1) = 0*L(ind1,ind1) - eta*1/dx;
    L(ind1,ind2) = 0*L(ind1,ind2) + eta*1/dx;
  end
end

% find the vector that is in the null space which we need to be able to
% compute the appropriate amount of flux out of the boundary
%tic
[u,s,v] = svds(L',1,'smallest');
%toc
%tic
%[V,D] = eig(full(L'));
%toc
%[~,ind] = min(abs(diag(D)));
%nullVec = V(:,ind);
nullVec = v;
%norm(nullVec2 - nullVec)
%norm(L'*nullVec)
%norm(L'*nullVec2)
%pause


end % PoissonSolverNeumannMatrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [psi,psi_x,psi_y] = PoissonSolverNeumann(o,f)
% Solve the 2D Poisson equation using centered differences and
% Neumann boundary conditions

N = size(f,1);
psi_x = zeros(N);
psi_y = zeros(N);

% remove the edge terms since we are assuming a homogeneous Neumann
% boundary condition
f = f(2:end-1,2:end-1); 

% choose the boundary data so that the compatability condition is
% satisfied. Checked with Fredholm Alternative and this worked exactly
% as we need (see lines below)
%bdFlux = +o.dx/(4*(N-2)) * sum(f(:));

rhs = zeros(N^2,1);
M = size(f,1);
for k = 1:M
  rhs(N+(k-1)*N+2:N+(k-1)*N+M+1) = f(:,k);
end

indLeft = o.strucLeft;
indBot = o.strucBot;
len = o.strucSize;
[cy,cx] = meshgrid(1:o.N,1:o.N);

for j = 1:numel(indLeft)
  % normal derivative at lower left corner
  ind = sub2ind([N,N],indLeft(j),indBot(j));
  rhs(ind) = +o.s*cos(o.windDir);

  % normal derivative at upper left corner
  ind = sub2ind([N,N],indLeft(j),indBot(j)+len);
  rhs(ind) = -o.s*cos(o.windDir);

  % normal derivative at lower right corner
  ind = sub2ind([N,N],indLeft(j)+len,indBot(j));
  rhs(ind) = +o.s*cos(o.windDir);

  % normal derivative at upper right corner
  ind = sub2ind([N,N],indLeft(j)+len,indBot(j)+len);
  rhs(ind) = -o.s*cos(o.windDir);

  % start of imposing Neumann boundary condition along the sides
  % (non-corners)
  for k = indLeft(j)+1:indLeft(j) + len - 1
    % update indices for normal derivative along bottom side
    ind = sub2ind([N,N],k,indBot(j));
    rhs(ind) = +o.s*cos(o.windDir);

    % update indices for normal derivative along top side
    ind = sub2ind([N,N],k,indBot(j)+len);
    rhs(ind) = -o.s*cos(o.windDir);
  end

  for k = indBot(j)+1:indBot(j) + len - 1
    % update indices for normal derivative along left side
    ind = sub2ind([N,N],indLeft(j),k);
    rhs(ind) = -o.s*sin(o.windDir);

    % update indices for normal derivative along right side
    ind = sub2ind([N,N],indLeft(j)+len,k);
    rhs(ind) = -o.s*sin(o.windDir);
  end
end


%[V,D] = eig(full(o.L)');
%[~,ind] = min(abs(diag(D)));
%
% geometric way of computing boundary flux, but not sure what to put in
% the case of obstacles
%bdFlux = sum(rhs)/(4*(o.N-2));

% linear algebraic way of computing boundary flux, but not sure how to
% come up with this number physically in the case of obstacles
%bdFlux = sum(rhs.*V(:,ind))/...
%    (sum(V(1:N,ind)) + sum(V(N^2-N+1:N^2,ind)) + ...
%    sum(V(N+1:N:N^2-2*N+1,ind)) + sum(V(2*N:N:N^2-N,ind)));
bdFlux = sum(rhs.*o.nullVec)/...
    (sum(o.nullVec(1:N)) + ...
     sum(o.nullVec(N^2-N+1:N^2)) + ...
     sum(o.nullVec(N+1:N:N^2-2*N+1)) + ...
     sum(o.nullVec(2*N:N:N^2-N)));


% put constant flux at all points on the four sides of the boundary
rhs(1:N) = -bdFlux;
rhs(N^2-N+1:N^2) = -bdFlux;
rhs(N+1:N:N^2-2*N+1) = -bdFlux;
rhs(2*N:N:N^2-N) = -bdFlux;

%sum(rhs.*V(:,ind))
%norm(o.L*(o.L\rhs) - rhs,inf)
%disp('here')
%pause

% solve for the potential
psi = o.L\rhs;

% reshape the potential so that it agrees with the mesh grid
psi = reshape(psi,N,N);

%% Use centered differencing to compute the gradient of the potential
%% (ie. the velocity)
%psi_x(2:end-1,2:end-1) = 0.5*(psi(3:end,2:end-1) - ...
%      psi(1:end-2,2:end-1))/o.dx;
%psi_y(2:end-1,2:end-1) = 0.5*(psi(2:end-1,3:end) - ...
%      psi(2:end-1,1:end-2))/o.dx;

[psi_y,psi_x] = gradient(psi);

%[cx(10:20,30) cy(10:20,30) psi(10:20,30) psi(10:20,30)-psi(10:20,29)]
%[cx(10:20,30) cy(10:20,30) psi(10:20,40) psi(10:20,41)-psi(10:20,40)]
%pause

% boundary conditions
psi_x(1,:) = -bdFlux;
psi_x(end,:) = +bdFlux;
psi_y(:,1) = -bdFlux;
psi_y(:,end) = +bdFlux;


end % PoissonSolverNeumann

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [psi,psi_x,psi_y] = PoissonSolverDirichlet(o,f)
% Solve the 2D Poisson equation using centered differences and
% homogeneous Dirichlet boundary conditions

N = size(f,1);
% remove the edge terms since we are assuming a Dirichlet boundary
% condition
f = f(2:end-1,2:end-1); 

e = ones(N-2,1);
A = spdiags([e -2*e e],(-1:1),N-2,N-2);
L = (kron(A,eye(N-2)) + kron(eye(N-2),A));

psi = zeros(N);
psi_x = zeros(N);
psi_y = zeros(N);

psi(2:end-1,2:end-1) = reshape(L\(o.dx^2*f(:)),N-2,N-2);
psi_x(2:end-1,2:end-1) = 0.5*(psi(3:end,2:end-1) - ...
      psi(1:end-2,2:end-1))/o.dx;
psi_y(2:end-1,2:end-1) = 0.5*(psi(2:end-1,3:end) - ...
      psi(2:end-1,1:end-2))/o.dx;

end % PoissonSolverDirichlet


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [velx,vely] = computeVelocity(o,psix,psiy,etax,etay)

velx = psix + etay + o.s*sin(o.windDir);  
vely = psiy - etax + o.s*cos(o.windDir); 

indLeft = o.strucLeft;
indBot = o.strucBot;
len = o.strucSize;

[cy,cx] = meshgrid(1:o.N,1:o.N);
%[cx(10:20,30) cy(10:20,30) psiy(10:20,31)]


% zero out interior points
for j = 1:numel(indLeft)
  velx(indLeft(j):indLeft(j)+len,indBot(j):indBot(j)+len) = 0;
  vely(indLeft(j):indLeft(j)+len,indBot(j):indBot(j)+len) = 0;
end
%clf
%quiver(cx,cy,velx,vely)
%pause

if o.emberWash
  % norm of the velocity
  normvel = sqrt(velx.^2 + vely.^2);
  % sine and cosine of the velocity field
  cosarg = velx./normvel;
  sinarg = vely./normvel;

  % mean of max distance travelled by an ember
  ExpProb = 2*o.s;
  % Probability that an ember is generated and creates a new fire
  BernoulliProb = 0.2;
  rdbino = binornd(1,BernoulliProb,o.N,o.N);
  % increase velocity with an exponential distribution in the direction
  % of the fire-induced and background wind velocity
  velx = velx + rdbino.*cosarg.*exprnd(ExpProb,o.N,o.N);
  vely = vely + rdbino.*sinarg.*exprnd(ExpProb,o.N,o.N);
end

end % computeVelocity

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function updateState(o,cells)

% index of points currently on fire
[FireIndx,FireIndy] = find(cells.state == 2);

% update burn time of combusting cells
for k = 1:numel(FireIndx)
  cells.burnTime(FireIndx(k),FireIndy(k)) = ...
        cells.burnTime(FireIndx(k),FireIndy(k)) + 1;
end

NewFireIndx = [];
NewFireIndy = [];
for k = 1:numel(FireIndx)
  ii = FireIndx(k); jj = FireIndy(k);
  % Compute index of points along the Bresenham line starting at burning
  % cells and going in direction of velocity vector
  [indx,indy] = cells.bham(cells.cx(ii,jj),cells.cy(ii,jj),...
        o.dt/o.dx*cells.velx(ii,jj),o.dt/o.dx*cells.vely(ii,jj));

  % first point is always the cell itself, so don't need to consider
  idx = indx(1); idy = indy(1); 
  indx = indx(2:end); % cells along bresenham line 
  indy = indy(2:end); % cells along bresenham line 

  % distance that an ember will travel is wind speed times a normal
  % random variable for the flight time. Also multiply by a bernoulli
  % since not all launched embers create a new fire 
  if o.spotting
    % mean and variance of flight time in units of seconds
    MeanFlightTime = 20;
    VarFlightTime = 5;
    % probability that an ember is generated and ignites a spot fire
    BernoulliProb = 0.01;
    EmberDist = o.s*normrnd(MeanFlightTime,VarFlightTime,1)*...
      binornd(1,BernoulliProb,1);
    if EmberDist > 0
      indx = [indx; idx + round(EmberDist*sin(o.windDir))];
      indy = [indy; idy + round(EmberDist*cos(o.windDir))];
    end
  end

  % Ignite along the Bresenham line but stop as soon as you reach a cell
  % that is either already burning or burnt out
  j = 1;
  while (j <= numel(indx) && ...
         indx(j) > 0 && indx(j) < cells.N && ...
         indy(j) > 0 && indy(j) < cells.N && ...
         cells.state(indx(j),indy(j)) == 1)
    % Save indicies of points that will be ignited
    NewFireIndx = [NewFireIndx indx(j)];
    NewFireIndy = [NewFireIndy indy(j)];
    j = j + 1;
  end

  if cells.burnTime(ii,jj) == cells.flameOut(ii,jj)
    % START OF DIFFUSION MODEL
    % diffuse cells to the east
    if (cells.state(ii+1,jj) == 1 && rand < o.probIgnite) 
      NewFireIndx = [NewFireIndx ii+1];
      NewFireIndy = [NewFireIndy jj+0];
    end
    % diffuse cells to the west
    if (cells.state(ii-1,jj) == 1 && rand < o.probIgnite) 
      NewFireIndx = [NewFireIndx ii-1];
      NewFireIndy = [NewFireIndy jj+0];
    end
    % diffuse cells to the north
    if (cells.state(ii,jj+1) == 1 && rand < o.probIgnite) 
      NewFireIndx = [NewFireIndx ii+0];
      NewFireIndy = [NewFireIndy jj+1];
    end
    % diffuse cells to the south
    if (cells.state(ii,jj-1) == 1 && rand < o.probIgnite) 
      NewFireIndx = [NewFireIndx ii+0];
      NewFireIndy = [NewFireIndy jj-1];
    end
    % diffuse cells to the northeast
    if (cells.state(ii+1,jj+1) == 1 && rand*sqrt(2) < o.probIgnite) 
      NewFireIndx = [NewFireIndx ii+1];
      NewFireIndy = [NewFireIndy jj+1];
    end
    % diffuse cells to the northwest
    if (cells.state(ii-1,jj+1) == 1 && rand*sqrt(2) < o.probIgnite) 
      NewFireIndx = [NewFireIndx ii-1];
      NewFireIndy = [NewFireIndy jj+1];
    end
    % diffuse cells to the southeast
    if (cells.state(ii+1,jj-1) == 1 && rand*sqrt(2) < o.probIgnite) 
      NewFireIndx = [NewFireIndx ii+1];
      NewFireIndy = [NewFireIndy jj-1];
    end
    % diffuse cells to the southwest
    if (cells.state(ii-1,jj-1) == 1 && rand*sqrt(2) < o.probIgnite) 
      NewFireIndx = [NewFireIndx ii-1];
      NewFireIndy = [NewFireIndy jj-1];
    end
    % END OF DIFFUSION MODEL
  end
end

% Ignite all new cells
for k = 1:numel(NewFireIndx)
  cells.state(NewFireIndx(k),NewFireIndy(k)) = 2;
end

% burn out cells that have been burning for the maximum number of
% time steps.
for k = 1:numel(FireIndx)
  if cells.burnTime(FireIndx(k),FireIndy(k)) == ...
        cells.flameOut(FireIndx(k),FireIndy(k));
    cells.state(FireIndx(k),FireIndy(k)) = 0;
  end
end



end % updateState

end % methods

end % classdef
