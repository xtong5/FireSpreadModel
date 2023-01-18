clear
for strength=[0, 25, 50, 100]
  folderName = "frames_" + int2str(strength);
  if ~exist(folderName, 'file')
    mkdir(folderName);
  end
  system('rm frames');
  system("ln -s " + folderName + " frames");
  fireLine
  clear
end