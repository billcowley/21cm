%  Two-freq spectrum calculations for H1 signals, using RPI/py frontend 
%  with rtlpars.txt file
% 
%
%  This script:  reads psdx.txt files as they become available 
%                reads text file called RA_settings.txt to get Az and El
%                store results in spect.log
%                

more off

if exist('p2res')~=1,  p2res = 1; end    % plot both the final estimates 

if exist('runcomment')~=1,
    runcomment='please insert notes here',
end


% get the key params from the parameter files...
pfname = 'rtlpars.txt';   pf = fopen(pfname, 'r');
if pf<0, error(' *** can not open parameter file '); end
pars = fscanf(pf,'%s');  fclose(pf);
eval(pars,'error(" *** param evaluation error")')


Vcol = ['r' 'y' 'm' 'g' 'c' 'b' 'k'];
Ncol = 7;   ncol = 1;   Mcol =12;

if (exist('PROCoctave')~=1),  PROCoctave=1,   end;   


if (exist('Nt')~=1)    Nt = 100000, end     % max loop count
if (exist('plt')~=1)   plt = 1,    end      % controls plots and debug

if (exist('logs')~=1)  logs=0,     end      % control spectra log file
if (exist('Nav')~=1)   Nav = 20;    end;    % number of spectral averages


if (exist('Ymin')~=1)  Ymin = -70,  end
if (exist('Yrng')~=1)  Yrng = 60,   end

if (exist('Tpause')~=1)    Tpause = 5,              end
if (exist('slogname')~=1)  slogname = 'spect.log',   end
if (exist('azelname')~=1)  azelname = 'RA_settings.txt',   end

if exist('Nskip')~=1, Nskip = Nw/4,     end;  % numb of bins to be skipped in

% following are read from pars file 
p.Fs = Fs;  p.Nw = Nw;
p.Fc1 = Fc1;   p.Fc2 = Fc2;
p.Gain = Gain;
frq_DC = (-Nw/2:Nw/2-1)*Fs/Nw;
logs=1;

if (logs==1)
    fslog = fopen(slogname, 'a');
    s= date;
    s2= datestr(now, 'HH:MM:SS');
    fprintf(fslog, 'Start spectra log on %s %s\n',s,s2);
    fprintf(fslog, 'params: %s \n',pars);
    fprintf(fslog, 'run notes: %s\n', runcomment)
    disp(['Logging spectra with params: ' pars]);

    fprintf(fslog, 'frequency bins: %s',s);
    fprintf(fslog, ' %.3f',(Fc1+frq_DC)/1e6);
    fprintf(fslog, '\n');
    fclose(fslog);
else
    disp('no spectra logging');
end

p.logs = logs;

nav = 0;

df2name = 'psd2.txt'
nloop=0;
if plt>0  figure(12); clf; end

while (1)
   af1 = stat('psd1.txt'); 
   af2 = stat(df2name);
   % 8 bytes per bin 
   if ((af1.size >=8*Nw) && (af2.size >=8*Nw))
        
        nav = nav + 1;
        pause(1);    % wait (?) 
        % now we assume both PSDs are ready to load in dBs
	newdata = 0;
	while (newdata==0)
  	  try
            %load '-ascii'  psd1.txt;
            %load '-ascii'  psd2.txt;
            p1 = fopen('psd1.txt'); 
            p2 = fopen('psd2.txt');
            clear ps1
	    clear ps2
            ps1 = textscan(p1, '%f');
            ps2 = textscan(p2, '%f');
            psd1 = cell2mat(ps1);
            psd2 = cell2mat(ps2);
            fclose all
	    if ((length(psd1) == Nw ) && (length(psd2)==Nw)) newdata = 1;
            else   disp('not enough psd values !');
            end
          catch
	     disp('load error from psd files');
             pause(10)
	     fclose all
	  end_try_catch
	end 
        psd1=fftshift(psd1); 
        psd2=fftshift(psd2);
        psdL1 = 10.^(psd1/10); psdL2 = 10.^(psd2/10); 
	
        % open Az/El file details
        fazel = fopen(azelname, 'r');
        if fazel<0,
            disp('*** could not open RA settings file, use dummies !! ');
            sal = 'please fix the RA settings !!!'
        else
            sal = fscanf(fazel, '%s'); 
            disp(['orientation is ' sal]); 
            fclose(fazel);
        end    
        if nav==1
            psd_sum1 = psdL1;   psd_sum2 = psdL2;
        else
            psd_sum1 = psd_sum1 + psdL1;   psd_sum2 = psd_sum2 + psdL2;
        end
        
	psdmax = max([psd_sum1' psd_sum2']); 
      
        psd_dB1 = 10*log10(psd_sum1/psdmax);
        psd_dB2 = 10*log10(psd_sum2/psdmax);
        psd_dB_diff = psd_dB1 - psd_dB2;
        s= datestr(now, 'HH:MM:SS');
        disp(['nav is ' num2str(nav) ' at ' s]);
        frq_range = (1+Nskip:length(frq_DC)-Nskip);

        if plt>0
            figure(9);  clf;
            plot((-1400e6+Fc1+frq_DC(1+Nskip: end-Nskip))/1e6, psd1(1+Nskip:end-Nskip));
            xlabel('Freq (MHz)');    ylabel('PSD (dB)');
            grid on;
            title(['First Spectral Estimate' ]);
            hold on;
          
            plot((-1400e6+Fc2+frq_DC(1+Nskip: end-Nskip))/1e6, psd2(1+Nskip:end-Nskip), 'g');
            xlabel('Freq (MHz)');    ylabel('RPI Spectra');
            grid on;
            title(['Both Spectral Estimates ']);
            
            
            figure(11);  clf;
            plot((-1400e6+Fc1+frq_DC(1+Nskip:end-Nskip))/1e6, psd_dB_diff(1+Nskip:end-Nskip));hold on;
            
            xlabel('Freq (MHz)');    ylabel('PSD (dB)');
            %s= datestr(now, 'HH:MM:SS');
            % disp(['nav is ' num2str(nav) ' at ' s]);
            title(['Spectral Difference #  ' num2str(nav) ' at ' s]);
            grid on;
        end
        
        if ((logs==1)&&(nav==Nav))
            disp('write results to log file');
            fslog = fopen(slogname, 'a');
            s= datestr(now, 'yyyymmddTHHMMSS');
            pmx = [' Pmx ' num2str(psdmax)];          % added 07042018 
            fprintf(fslog, '%s %s %s',s, sal, pmx);
            fprintf(fslog, '\n');
            fprintf(fslog, ' %.2f',psd_dB1);
            fprintf(fslog, '\n');
            fprintf(fslog, ' %.2f',psd_dB2);
            fprintf(fslog, '\n')
            fclose(fslog);
	end
            
            %fflog = fopen(flogname, 'a');    % final results file
            %fprintf(fflog, '%s ',sal);
            %fprintf(fflog, '%s,%s,',s,sal);
            %fprintf(fflog, ' %.2f',psd_dB_diff);
            %fprintf(fflog, '\n');
            %fclose(fflog);      

            % frq_range = (1+Nskip:length(frq_DC)-Nskip);

         if ((nav==Nav) && (plt>0)) 
  	       ncol=ncol+1;
	       if ncol>Ncol, ncol=1; end
               figure(12);
               plot((Fc1+frq_DC(frq_range))/1e6, ...
                psd_dB_diff(frq_range),  Vcol(ncol));
               xlabel('Freq (MHz)');    ylabel('PSD (dB)');
               s= datestr(now, 'HH:MM:SS');   % disp(s);
               title(['Final Spectral Difference: + version at ' s]);
            
               if p2res
                hold on; 
                plot((Fc2+frq_DC(frq_range))/1e6, ...
                    -psd_dB_diff(frq_range),  Vcol(ncol));
                xlabel('Freq (MHz)');    ylabel('PSD (dB)');
                s= datestr(now, 'HH:MM:SS');   % disp(s);
                title(['Final Spectral Difference: +/- version at ' s]);
               end
	       hold on 
	  end
            
           
        
        if nav==Nav,
	  nav = 0; nloop=nloop+1;
          clear psdL1 psdL2
          [maxX, maxI] = max(psd_dB_diff);
          maxFreq = Fc1+ frq_DC(maxI);
          disp(['  ****  Location of largest peak is ' num2str(maxFreq/1e6) ' MHz'])
          pk2pk = (max(psd_dB_diff(frq_range))- min(psd_dB_diff(frq_range)))/2;
          disp(['  >>>  Pk max = ' num2str(pk2pk,2) ' dB']);
        end
        
        fflush(stdout); pause(0.1);
        ch = kbhit(1);   % 
        if (length(ch)>0)  break;   end
    
       
    else
        disp( 'waiting for psd2 file from python')
        pause(1)
    end
    
    disp([' loop '   num2str(nloop)]);
    
    disp(['pausing for ' num2str(Tpause) ' secs']);
    pause(Tpause);
    
end
