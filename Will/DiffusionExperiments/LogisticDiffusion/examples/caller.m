clear
for strength=[0, 0.01, 0.05, 0.1, 0.25, 0.5, 1, 5, 10, 50, 100, 500, 1000]
  folderName = sprintf('frames_%.2f',strength);
  if ~exist(folderName, 'file')
    mkdir(folderName);
  end
  system('rm frames');
  system("ln -s " + folderName + " frames");
  fireLine
  clear
end

