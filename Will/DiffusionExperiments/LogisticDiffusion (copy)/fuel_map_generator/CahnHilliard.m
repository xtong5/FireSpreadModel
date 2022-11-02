Ngrid = 401;
N = 2^ceil(log2(Ngrid));

x = (0:N-1)';
[x,y] = meshgrid(x,x);

u = zeros(size(x));
prop = 0.5;
p = randperm(N^2);
u(p(1:ceil(prop*N^2))) = 1;
u(p(1+ceil(prop*N^2):end)) = -1;
uold = u;

%u = 2*rand(N)-1 + 0.2;
surf(x,y,u)
view(2); shading interp;
axis equal
axis([0 Ngrid 0 Ngrid])
pause(.01)
%pause
s = sum(sum(u));

modes = (-(N-1)/2:(N-1)/2)';
[modesx,modesy] = meshgrid(modes,modes);
err = [];

dt = 1e-3;
eps = 1e-2; % parameter that determines the correlation length scale
for k = 1:1000
  fu = u.^3 - u;
  fuh = fftshift(fft2(fu));
  Lfuh = -(modesx.^2 + modesy.^2).*fuh;
  rhsh = fftshift(fft2(u)) + dt*Lfuh;
  uh = rhsh./(1 + eps*dt*(modesx.^4 + modesy.^4)); 
  u = real(ifft2(ifftshift(uh)));

  s = [s sum(sum(u))];
  err = [err norm(u(:) - uold(:),inf)];
  uold = u;

  if 1
  subplot(2,1,1)
  surf(x,y,u);colorbar; caxis([-1 1])
  view(2); shading interp;
  axis equal
  axis([0 Ngrid 0 Ngrid])
  subplot(2,1,2);
  semilogy(err)
  pause(.01)
  end
end


