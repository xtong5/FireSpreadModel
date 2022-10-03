clear
for strength=[1, 10, 100]
  system('rm frames')
  system(strcat('ln -s frames_', int2str(strength), ' frames'))
  fireLine
  clear
end