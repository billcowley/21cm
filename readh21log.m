
% read and display the H21 data using log format from 2018, Bill
%
% the mesh plot can be rotated using the figure viewing tools
% this version tries to handle multiple sets of data in the one log file!
%
% set Npt number of spectra per plot;
% or set numb of hours to plots ('samePH')
% add 2nd plot showing results in 1D
% use deb = 2 or 3 to plot raw spectral for each record
%
% add option for movie of the raw PSDs -- only on matlab version !

if exist('deb')~=1,    deb =0,    end
if exist('Npt')~=1,    Npt=2000,   end;  % how many records per plot?
if exist('plt')~=1,    plt =1,    end;  % one plot or 2 plots?
if exist('Fmin')~=1,   Fmin=1420.1; end  % limit the plot range
if exist('Fmax')~=1,   Fmax=1420.8; end
if exist('samePH')~=1, samePH =6, end % plot period in hrs
if exist('Iskip')~=1,  Iskip = 2; end % initial records to skip in plot
if exist('psdVid')~=1, psdVid =0; end % set to 1 for video of raw PSDs
if exist('Nvid')~=1,   Nvid = 1500; end;  % frames in video
if exist('poly7')~=1,  poly7 =0;    end   % used for poly fitting method
                                          % TBC  see code below at EOF
if exist('chkGIF')~=1, chkGIF=0;    end   % use for making animated gif 

samesecs = samePH * 3600;      % conv to secs

% if exist('noslope')~=1,  noslope = 0; slop=[]; end % controls PSD slope adj
% if noslope
%     if exist('Fit1')~=1,  Fit1= 1420.1e6, end
%     if exist('Fit2')~=1,  Fit2= 1420.9e6, end
%     if exist('Fitw')~=1,  Fitw= 0.15e6,   end
% end

if exist('fname')==0, fname = 'spect.log',  end
%eg  fname = 'spect22_24.log'
fid = fopen(fname);

if (fid==0) error([' *** could not open log file: ' fname]); end
disp(['opened log file for H21 project: ' fname]);



if psdVid==1,
    v=version;
    if v(1)=='4',  % assume octave
        disp('no movie with octave ! ');
        psdVid=0;
    else
        jvid = 0;
        vidObj = VideoWriter('h21psd.avi');
        disp('created video object ');
        open(vidObj);
    end
end

while (1)
    n =1; clear pows; clear secs pmax;  clear pks*; sameperiod =1;
    
    while (~feof(fid))&&(n<Npt) % while not EOF
        
        line1 = fgets(fid);
        
        if line1(1:5)=='Start'
            disp(line1);
        elseif line1(1:7)=='params:'
            disp(line1);
            line1(1:7)=[];    % delete initial chars
            eval(line1,'error(" *** param evaluation error")')
            frq = (Fc1+ (-Nw/2:Nw/2-1)*Fs/Nw)/ 1e6;
            if (poly7==1)
                polyfile = 'poly_h21_512.mat'
                load(polyfile);
                if (fit.Fc1==Fc1  && fit.Fc2==Fc2 && fit.Fs==Fs && fit.Nw==Nw)
                    disp('loaded poly fits ');
                else
                    disp('wrong param values fit model- revert to subtraction')
                    poly7=0,
                end
            end
            
        elseif line1(9)=='T'
            if deb>1, disp(line1); end
            secs(n) = 3600*str2num(line1(10:11)) + ...
                60*str2num(line1(12:13)) + str2num(line1(14:15));
            if n>1
                if floor(secs(n)/samesecs) > floor(secs(n-1)/samesecs)
                    sameperiod =0;
                end
            end
            timenow = line1(1:15);
            if n==1,
                startup = line1(1:15);
                startup(9) ='-';
                disp(['starting at ' startup]);
            end;
            inx=findstr(line1, 'Pmx');
            if (inx>0)
                pmax(n) = str2num(line1(inx+3:end));
            end
            
        elseif line1(1:4)=='freq'
            disp(['start of frequency bins list: ' line1(1:24)]);
            line1(1:24)=[];
            frq = sscanf(line1, ' %f', Nw);
        else
            pow1= sscanf(line1, ' %f', Nw);
            line2 = fgets(fid);
            pow2 = sscanf(line2, ' %f', Nw);
            if ((length(pow1)==Nw) & (length(pow2)==Nw))
                if poly7==0
                    psd_diff = pow1'-pow2';    % store pow in array
                else
                    pp1 = pow1' -polyval(pol, frq-1420);
                    pp2 = pow2' - polyval(pol2, frq-1420);
                    k_polm = round((Fc2-Fc1)/(Fs/Nw));
                    pp3 = pp1(k_polm:Nw) + pp2(1:Nw-k_polm+1);
                    psd_diff = pp3;
                    
                end
                pows(n, :)= psd_diff;
                [maxx, maxi] = max(pows(n,:));
                pks_max(n) = maxx;   pks_frq(n) = frq(maxi);
                if deb>1,
                    disp([' read two spectra for t = ' num2str(secs(n))])
                end
                
                if deb>0
                    if deb==2
                        f201=figure(201); clf
                        plot(pow1); hold on;
                        plot(pow2, 'k');
                        title(timenow);
                        ax = axis; ax(3) = ax(4)-5; % 5 dB range
                        axis(ax);
                        pause(0.02);
                        
                        if psdVid==1
                            jvid = jvid + 1;
                            currFrame= getframe(201);
                            writeVideo(vidObj,currFrame);
                            if jvid > Nvid
                                psdVid =0;
                                close(vidObj);
                                disp('closed movie')
                            end
                        end
                    end
                    
                    if (deb==6) pause; end
                    if deb==4
                        f202=figure(202); clf
                        plot(frq, pow1); hold on;
                        plot(frq+(Fc2-Fc1)/1e6, pow2, 'k');
                        title(timenow);
                        ax = axis; ax(3) = ax(4)-5; % 5 dB range
                        axis(ax);
                        pause(0.02);
                    end
                    if deb==5
                        f203=figure(203); clf
                        plot(frq/1e3, psd_diff); hold on;
                        
                        title(timenow);
                        ax = axis; ax(3) = ax(4)-5; % 5 dB range
                        axis(ax);
                        pause(0.02);
                    end
                    if chkGIF>0,
                        h_now = str2num(timenow(end-5:end-4)); % get hours now
                        if ((makeGIF.start>=h_now) & (makeGIF.stop>h_now))
                            if makeGIF.active==0,
                                makeGIF.active =1;  
                                if chkGIF==201, gif('h21_anim.gif', 'frame', f201);deb=2; end
                                if chkGIF==202, gif('h21_anim.gif', 'frame', f202);deb=4; end
                                if chkGIF==203, gif('h21_anim.gif', 'frame', f203);deb=5; end
                            else
                                if chkGIF==201, gif('frame', f201); end
                                if chkGIF==202, gif('frame', f202); end
                                if chkGIF==203, gif('frame', f203); end
                            end
                        else
                            if makeGIF.active ==1,
                                makeGIF.active = 0; 
                                %gif('clear', 'frame', f201); 
                                deb = 1;
                            end
                        end
                    end
                    
                end
                n=n+1;
                if (sameperiod==0) break; end;
            else
                disp('pow1 and pow2 are different lengths??')
            end
        end
        
    end
    
    
    if length(secs)>1
        disp(['read ' num2str(length(secs)) ' records']);
        sec0=0; days =0;
        secn = zeros(1,length(secs));
        for i =1: length(secs)
            if secs(i)<sec0,
                disp('next day found')
                days=days+1;
            end
            sec0 = secs(i);
            secn(i) = secs(i)+days*24*3600;
        end
        
        %         nxti = secs<10000;   % fix after midnight times
        %         secs(nxti) = secs(nxti) + 24*3600;
        
        %make mesh plot
        figure(100)
        disp('plotting results in mesh plot ...');
        if poly7==1,
            frqp = frq(k_polm:Nw); pows=pows(:, 1:Nw-k_polm+1);
            f_ind = frqp>Fmin & frqp<Fmax;
            frq2 = frqp(f_ind);
        else
            f_ind = frq>Fmin & frq<Fmax;
            frq2 = frq(f_ind);
        end
        
        pows=pows(:, f_ind);
        sz = size(pows);
        
        mesh(frq2, secn(Iskip:sz(1))/3600, pows(Iskip:end, :))
        title(['H21 Spectral Estimates from log file ' fname], 'Interpreter','none')
        xlabel('Frequency in MHz')
        zlabel('Power (dB)')
        ylab = ['Hours; started at ' startup];
        ylabel(ylab)
        
        if exist('pmax')==1
            figure(105);
            plot(secs/3600, 10*log10(pmax));
            xlabel('Hours'); ylabel('Max Power')
        end
        
        if plt>1
            % show intensity and colour in 1D
            disp('replot in 1D')
            figure(101)
            %cmap=colormap('default');
            cmap = colormap;
            Am= min(pks_max); AM=max(pks_max);
            %Fm= min(pks_frq); FM=max(pks_frq);
            Fm = 1420.2; FM = 1420.8;
            for n=1:length(secn)
                h=plot(secn(n),0, 'x'); hold on;
                nw = 1+floor((pks_max(n)-Am)/(AM-Am)*15);
                nc =floor((pks_frq(n)-Fm)/(FM-Fm)*63);
                if nc<1, nc=1; end;  if nc>63, nc=63; end
                set(h,'color', cmap(nc,:), 'markersize',nw);
                title('1D peak plot: size=level, colour=freq')
                xlabel(ylab);
            end
        end
        disp('pausing until key pressed');
        if (psdVid==0)  pause;
        else  pause(3);  end
    end
    if feof(fid), disp('all done'); break, end
    
end
fclose(fid);



%{
% commands used to model no-signal psd for 'pow1'  at 1420MHz 
% saved with 
% --- save('poly_h21_512.mat', 'pol', 'pol2', 'fit', '-v7');

rng = 101:400;
figure;   plot(frq(rng), pow1(rng))

fpp = frq(rng);    ppp = pow1(rng);


% following not needed -- DC offset has been removed 
%ppp(158)=-.1 
%ppp(156)=-.05
%ppp(157)=-.09

pol = polyfit(frq(rng)-1420, pow1(rng)', 9);

pmod = polyval(pol, fpp-1420);
figure; plot(fpp, pmod)

plot(frq(rng), pmod, 'g')
%}
