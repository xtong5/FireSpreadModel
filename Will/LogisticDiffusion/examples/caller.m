% calls fireLine.m for multiple values of a parameter
% which must be handled correctly in fireLine
clear
for strength=[100, 150, 200, 250, 300, 350, 400, 450, 500]
  folderName = sprintf('frames_%.2f',strength);
  if ~exist(folderName, 'file')
    mkdir(folderName);
  end
  system('rm frames');
  system("ln -s " + folderName + " frames");
  fireLine
  clear
end
