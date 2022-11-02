function map = randTerrain(numAnchors, N)
  rng(1)
  anchors = rand(numAnchors, numAnchors);
  [x,y] = meshgrid(1:numAnchors);
  [xq,yq] = meshgrid(linspace(1,numAnchors,N));
  map = interp2(x,y,anchors,xq,yq,'cubic');
  map = map / 2;
  map = map + 0.5;
end