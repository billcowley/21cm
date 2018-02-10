% test code to rtl_proc.py :   just grab PSDs files and display
% bill feb 2018 

more off
pfname = 'rtlpars.txt';   pf = fopen(pfname, 'r');
if pf<0, error(' *** can not open parameter file '); end
pars = fscanf(pf,'%s');  fclose(pf);
disp(['params are ' pars])
eval(pars,'error(" *** param evaluation error")');

p.Fs = Fs;  p.Nw = Nw;   p.Nw = Nw; 
p.Fc1 = Fc1;   p.Fc2 = Fc2;
p.Gain = Gain;

frq_DC = fftshift(-Nw/2:Nw/2-1)*Fs/Nw; 

while (1)
  if length(stat('psd1.txt'))>0
  pause(0.1);
  load psd1.txt
  system('mv psd1.txt  psd1.old');
  figure(51);
  clf;
  plot(frq_DC+Fc1, psd1, 'b'); hold on
  printf("1")
  end
  
  if length(stat('psd2.txt'))>0
  pause(0.1);
  load psd2.txt
  system('mv psd2.txt  psd2.old');
  figure(51);
  plot(frq_DC+Fc2, psd2, 'k');
  printf("2")
  end
  
  pause(1)
  printf(".")

end  
