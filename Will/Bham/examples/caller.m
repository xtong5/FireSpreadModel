% calls fireLine.m for multiple values of a parameter
% which must be handled correctly in fireLine
clear
for strength=[0, 1000, 1100, 1200, 1300, 1400, 1500, 2000, 3000, 4000]
  frameFolder = sprintf('frames_%.2f',strength);
  if ~exist(frameFolder, 'file')
    mkdir(frameFolder);
  end
  system('rm frames');
  system("ln -s " + frameFolder + " frames");
  rng(0);
  fireLine
  dataFolder = sprintf('data_%.2f',strength);
%   if ~exist("dataFolder", 'file')
%     mkdir(dataFolder);
%     system(strcat("mv *.mat ", dataFolder));
%   end
  clear
end
