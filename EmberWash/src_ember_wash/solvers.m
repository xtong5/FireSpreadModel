classdef solvers
% different PDE solvers


properties
N
dx
dt
uwind
s
sinkStrength
vortStrength
probIgnite
L

end % properties

methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = solvers(prams)

o.N = prams.N;
o.dx = prams.dx;
o.dt = prams.dt;
o.uwind = prams.uwind;
o.s = prams.s;
o.sinkStrength = prams.sinkStrength;
o.vortStrength = prams.vortStrength;
o.probIgnite = prams.probIgnite;
o.L = o.PoissonSolverNeumannMatrix;


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
function L = PoissonSolverNeumannMatrix(o)
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
% Force the first term in psi to be 0
%L = [L;speye(1,N^2)];

% do LU decomposition so inversion can done efficiently
%[L,U] = lu(L);

end % PoissonSolverNeumannMatrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [psi,psi_x,psi_y] = PoissonSolverNeumann(o,f)
% Solve the 2D Poisson equation using centered differences and
% homogeneous Dirichlet boundary conditions

N = size(f,1);
psi_x = zeros(N);
psi_y = zeros(N);

% remove the edge terms since we are assuming a homogeneous Neumann
% boundary condition
f = f(2:end-1,2:end-1); 

% choose the boundary data so that the compatability condition is
% satisfied. Checked with Fredholm Alternative and this worked exactly
% as we need (see lines below)
bdFlux = +o.dx/(4*(N-2)) * sum(f(:));

rhs = zeros(N^2,1);
% put constant flux at all points on boundary
rhs(1:N) = -bdFlux;
rhs(N^2-N+1:N^2) = -bdFlux;
rhs(N+1:N:N^2-2*N+1) = -bdFlux;
rhs(2*N:N:N^2-N) = -bdFlux;
M = size(f,1);
for k = 1:M
  rhs(N+(k-1)*N+2:N+(k-1)*N+M+1) = f(:,k);
end
%rhs = [rhs;0];
psi = o.L\rhs;

%clf
%[V,D] = eig(full(o.L)');
%[~,s] = min(abs(diag(D)));
%V(:,s)'*rhs % by Fredholm Alternative, this number should be 0
%plot(o.L*psi - rhs)
%pause

psi = reshape(psi,N,N);

psi_x(2:end-1,2:end-1) = 0.5*(psi(3:end,2:end-1) - ...
      psi(1:end-2,2:end-1))/o.dx;
psi_y(2:end-1,2:end-1) = 0.5*(psi(2:end-1,3:end) - ...
      psi(2:end-1,1:end-2))/o.dx;

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

end % PoissonSolver


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [velx,vely] = computeVelocity(o,psix,psiy,etax,etay,cos,sin)

% rdbino = binornd(1,0.01,o.N,o.N);
rdbino = zeros(o.N,o.N);
velx = psix + etay + cos.*exprnd(o.s,o.N,o.N).*rdbino;
vely = psiy - etax + o.s + sin.*exprnd(o.s,o.N,o.N).*rdbino;
% velx = psix + etay + exprnd(10*o.s,o.N,o.N);
% vely = psiy - etax + o.s + exprnd(10*o.s,o.N,o.N); %v*arrival time 

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

%  [indx-0.0 indy-0.0]
%  pause

  % first point is always the cell itself, so don't need to consider
  indx = indx(2:end); % cells along bresenham line 
  indy = indy(2:end); % cells along bresenham line 

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

%  for j = 1:numel(indx)
%    % Save indicies of points that will be ignited
%    if (indx(j) > 0 && indx(j) < cells.N && ...
%        indy(j) > 0 && indy(j) < cells.N && ...
%        cells.state(indx(j),indy(j)) == 1)
%      NewFireIndx = [NewFireIndx indx(j)];
%      NewFireIndy = [NewFireIndy indy(j)];
%%      [ii jj indx(j) indy(j)]
%%      pause
%    end
%  end

  if cells.burnTime(ii,jj) == cells.flameOut
    % START OF DIFFUSION MODEL
    % diffuse cells to the east
    if (cells.state(ii+1,jj) == 1 && rand < o.probIgnite) 
  %    disp('east')
      NewFireIndx = [NewFireIndx ii+1];
      NewFireIndy = [NewFireIndy jj+0];
    end
    % diffuse cells to the west
    if (cells.state(ii-1,jj) == 1 && rand < o.probIgnite) 
  %    disp('west')
      NewFireIndx = [NewFireIndx ii-1];
      NewFireIndy = [NewFireIndy jj+0];
    end
    % diffuse cells to the north
    if (cells.state(ii,jj+1) == 1 && rand < o.probIgnite) 
  %    disp('north')
      NewFireIndx = [NewFireIndx ii+0];
      NewFireIndy = [NewFireIndy jj+1];
    end
    % diffuse cells to the south
    if (cells.state(ii,jj-1) == 1 && rand < o.probIgnite) 
  %    disp('south')
      NewFireIndx = [NewFireIndx ii+0];
      NewFireIndy = [NewFireIndy jj-1];
    end
    % diffuse cells to the northeast
    if (cells.state(ii+1,jj+1) == 1 && rand*sqrt(2) < o.probIgnite) 
  %    disp('northeast')
      NewFireIndx = [NewFireIndx ii+1];
      NewFireIndy = [NewFireIndy jj+1];
    end
    % diffuse cells to the northwest
    if (cells.state(ii-1,jj+1) == 1 && rand*sqrt(2) < o.probIgnite) 
  %    disp('northwest')
      NewFireIndx = [NewFireIndx ii-1];
      NewFireIndy = [NewFireIndy jj+1];
    end
    % diffuse cells to the southeast
    if (cells.state(ii+1,jj-1) == 1 && rand*sqrt(2) < o.probIgnite) 
  %    disp('southeast')
      NewFireIndx = [NewFireIndx ii+1];
      NewFireIndy = [NewFireIndy jj-1];
    end
    % diffuse cells to the southwest
    if (cells.state(ii-1,jj-1) == 1 && rand*sqrt(2) < o.probIgnite) 
  %    disp('southwest')
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
cells.state(cells.burnTime == cells.flameOut) = 0;

end % updateState

end % methods

end % classdef
