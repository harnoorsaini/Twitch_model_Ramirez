% =========================================================================
% Activation "twitch" Model as implemented in Muscle VUMAT 
% =========================================================================
% Harnoor Saini 
% Sep-2017
%
% Continuum Biomechanics and Mechanobiology 
% University of Stuttgart
% Stuttgart, Germany
%
% -------------------------------------------------------------------------
% convolutoion algorithm based on https://gist.github.com ... 
% /jshahbazi/5289503#file-fortran_convolution-f90
% -------------------------------------------------------------------------
% TO DOS/COMMENTS
% - f_sat does not work! (may leave it for now)
% - effect of frequency on amplitude is really strong - is this realistic?
% -------------------------------------------------------------------------
% RESET WORKSPACE
clear 
%close all
clc

outStr = ['INTEGRATED ACTIVATION BASED ON RAMIREZ 2012.'];
disp(outStr)
outStr = ['-------------------------------------------------------------'];
disp(outStr)

% -------------------------------------------------------------------------
% USER INPUTS & PARAMETERS
twocolumns = 1;
fatigue = 0;

% MU firing parameters
firing_type = 'generate'; % read or generate
firing_freq = 10; %Hz
filter_freq = 10; %Hz
firing_time = 0.5; %s

% unit pulse parameters 
Pprime = 0.11;
pulse_inc = 200;
pulse_tstep = 0.002; %s
Tcprime = 0.04; %s

% -------------------------------------------------------------------------
% SET UP FIRING TIMES & INITIALISE ARRAYS
% set up convolution or "total" (minimum) time required
firing_inc = firing_time/pulse_tstep;
conv_inc = firing_inc+pulse_inc-1;
conv_tvec = 0:pulse_tstep:pulse_tstep*(conv_inc-1);

% motor unit firings
switch firing_type
    case 'read'
        fname = ['C:\Users\saini\Documents\PhD_Local_C\41_MATLAB'...
            '\muscle_modelling\Twitch_model_Ramirez\code\mufirings.txt'];
        fileID = fopen(fname,'r');
        formatSpec = '%f';
        xMusf = fscanf(fileID,formatSpec);
        xMusf = xMusf';
        firing_inc = size(xMusf,2);
        conv_inc = firing_inc+pulse_inc-1;
        conv_tvec = 0:pulse_tstep:pulse_tstep*(conv_inc-1);
        fclose(fileID);
    case 'generate'
        firing_tstep = 1/firing_freq;
        for global_time = conv_tvec
            % if global_time is grea
            if global_time >= firing_tstep 
                t_offset = global_time - firing_tstep;
                break 
            end
        end
        if t_offset < pulse_tstep/2
            firing_tstep = firing_tstep - t_offset;
        elseif t_offset > pulse_tstep/2
            firing_tstep = global_time;
        else
            error('WARNING: firing_tstep = 0.5*pulse_tstep')
        end
        xMusf_temp = mod(int16(conv_tvec/pulse_tstep), ... 
            int16(firing_tstep/pulse_tstep))==0;
        xMusf(1:int16(firing_inc)) = xMusf_temp(1:int16(firing_inc));  
        if firing_tstep < pulse_tstep
            outStr = ['WARNING: increase pulse_tstep >=' num2str(firing_tstep)];
            disp(outStr)
        end

        % find the resulting effective firing frequency
        for i = 2:conv_inc
            if xMusf_temp(i) == 1
                eff_firing_freq = 1/conv_tvec(i);
                break
            end
        end

        if (eff_firing_freq ~= firing_freq) 
            outStr = ['WARNING: effective frequency is ' num2str(eff_firing_freq)];
            disp(outStr)   
            outStr = ['Desired frequency is ' num2str(firing_freq) ];
            disp(outStr)
            outStr = ['Consider reducing pulse_tstep'];
            disp(outStr)
        end  
end

% -------------------------------------------------------------------------
% COMPUTE THE UNIT PULSE
t = 0;
%Pprime = 5.0;
f_pulse(pulse_inc) = zeros;


for j = 1:pulse_inc
    f_pulse(j) = Pprime*t/Tcprime * exp(1-t/Tcprime); 
    t = t + pulse_tstep;
end

% -------------------------------------------------------------------------
% FREQUENCY DEPENDECE IN FORCE LEVEL
fr = firing_freq;
rfr = 1.0535;
cfr = 1.1245;
F_rnorm = fr*Tcprime;
F_fr = 1 - rfr * exp(-F_rnorm)/cfr;

if fatigue
    h = 3.5;
    A_m = 5.0;
    A_max = 0.5;
    A_0 = 0.0;
    t_Amax = 0.21;
    for j = 1:conv_inc
        A(j) = A_0 + (A_max-A_0) * (t/t_Amax) * exp((1-t)/t_Amax);
        f_sat(j) = (A(j)^h)/(A_m^h+A(j)^h); 
        A_max = A_max - 0.0002;
    end
else 
    f_sat = 1;
end

% -------------------------------------------------------------------------
% CONVOLUTE THE UNIT PULSE WITH THE MU FIRINGS
x = f_pulse; 
h = xMusf;

m = length(x);
n = length(h);

X = [x,zeros(1,n)]; 
H = [h,zeros(1,m)]; 

for i=1:n+m-1
    Y(i)=0;
    for j=1:m
        if(i-j+1>0)
                Y(i)=Y(i)+X(j)*H(i-j+1);
        else
        end
    end
end



% pad firings with 0s (for plotting)
for j = firing_inc+1:conv_inc
      xMusf(j) = 0;
end 

% global activations
alpha = Y;
alpha_ff = Y.*F_fr;
if fatigue, alpha_ff_fsat = Y.*f_sat*F_fr; end
% compare to inbuilt convolution
alpha_inbuilt = conv(f_pulse,double(xMusf));

d = designfilt('lowpassiir','FilterOrder',12, ...
    'HalfPowerFrequency',filter_freq/fr,'DesignMethod','butter');
alpha_filt = filtfilt(d,alpha);

figure(1)

hold on
plot(conv_tvec,alpha)
plot(conv_tvec,alpha_filt,'r')
plot(conv_tvec,alpha_ff,'g')
if fatigue, plot(conv_tvec,alpha_ff_fsat,'b'), end

%ylim([0 2])
plot(conv_tvec,alpha_inbuilt(1:length(conv_tvec)),'-b')
ylabel('Summed global activtation')
xlabel('time (s)')
%yyaxis right
%ylim([0 2])
%xlim([0 1])
stem(conv_tvec,xMusf)
%ylabel('MUAP firing')
title('MUAP and Summed Global Activation')


if fatigue, figure(2), plot(f_sat), end


fileID = fopen('2.txt','w');
outMat(1,:) = conv_tvec;
for i = 1:length(alpha_filt)
    if alpha_filt(i) < 0
        alpha_filt(i) = 0;
    end
end
outMat(2,:) = alpha_filt ;
if twocolumns, fprintf(fileID,'%12.8f %12.8f\n', outMat); end
if ~twocolumns
    for i = 1:length(alpha_filt)
        outStr = ['alpha_r(' num2str(i) ',1)=' num2str(conv_tvec(i))];
        fprintf(fileID,'%s\n', outStr);
    end
    fprintf(fileID,'%s\n', '-----');
    for i = 1:length(alpha_filt)
        outStr = ['alpha_r(' num2str(i) ',2)=' num2str(alpha_filt(i))];
        fprintf(fileID,'%s\n', outStr);    
    end
end
fclose(fileID);

disp('COMPLETE')

%hold on
%figure(2)
%pulse_tvec = 0:pulse_tstep:pulse_tstep*(pulse_inc-1);
%plot(pulse_tvec,f_pulse)

%figure(3)
%plot(0:1:100,F_fr)
%hold on

%figure(4)
%plot(conv_tvec,f_sat)
%hold on
%plot(conv_tvec,A)
