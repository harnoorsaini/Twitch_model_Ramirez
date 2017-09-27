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

% -------------------------------------------------------------------------
% USER INPUTS & PARAMETERS
% MU firing parameters
firing_type = 'generate'; % read or generate
firing_inc = 100; 
firing_freq = 20; %Hz (note 50Hz is tetanus)

% unit pulse parameters 
Pprime = 0.11;
pulse_inc = 200;
r = 1E7;
s = 1;
pulse_tstep = 0.005; %r*1E-9; %s
Tcprime = 0.04; %s*r*4E-9; %s

% set up convolution or "total" (minimum) time required
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
        firing_inc = size(xMusf,1);
        fclose(fileID);
    case 'generate'
        firing_tstep = 1/firing_freq;
        for global_time = conv_tvec
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
            error('warning, firing_tstep = 0.5*pulse_tstep')
        end
        xMusf_temp = mod(int8(conv_tvec/pulse_tstep), ... 
            int8(firing_tstep/pulse_tstep))==0;
        xMusf(1:firing_inc) = xMusf_temp(1:firing_inc);  
end

if firing_tstep < pulse_tstep
    outStr = ['WARNING: increase pulse_tstep >=' num2str(firing_tstep)];
    disp(outStr)

end


% -------------------------------------------------------------------------
% COMPUTE THE UNIT PULSE
t = 0;
%Pprime = 5.0;
f_pulse(pulse_inc) = zeros;

h = 1.05;
A_m = 0.55;
A_max = 0.5;
A_0 = 0;
t_Amax = 0.040;


for j = 1:pulse_inc
    A(j) = A_0 + (A_max-A_0) * (t/t_Amax) * exp((1-t)/t_Amax);
    f_sat(j) = A(j)^h/(A_m^h+A(j)^h); 
    A_max = A_max - 0.002;
    f_pulse(j) = Pprime*t/Tcprime * exp(1-t/Tcprime); 
    t = t + pulse_tstep;
end 

% -------------------------------------------------------------------------
% CONVOLUTE THE UNIT PULSE WITH THE MU FIRINGS
% pad firings with 0s
for j = firing_inc+1:conv_inc
      xMusf(j) = 0;
end 


% size of xw = length(f_pulse)+length(MF_firing)-1
m = 1;
w_max = 0.0;
summedFirings(conv_inc) = zeros;
for k = 1:conv_inc
      summedFirings(m) = 0.0;
% size of f_pulse
      for j = 1:firing_inc
            if ( (m-j+1) >= 1 )
                  summedFirings(m) = summedFirings(m) ... 
                      + f_pulse(j) * xMusf(m-j+1) * f_sat(j);
            end
      end 
      if (summedFirings(m) > w_max)
            w_max = summedFirings(m);
      end
      m = m + 1;  
end

% normalise activation
%alpha = summedFirings/w_max;
alpha = summedFirings;

% -------------------------------------------------------------------------
% FREQUENCY DEPENDECE IN FORCE LEVEL
rfr = 1.0535;
cfr = 1.1245;
F_rnorm = (0:1:100)*Tcprime;
F_fr = 1 - rfr * exp(-F_rnorm)/cfr;


% -------------------------------------------------------------------------
% SATURATION IN FORCE LEVEL
%h = 1.05;
%A_m = 0.55;
%A_max = 0.5;
%A_0 = 0;
%t_Amax = 0.040;


%t = 0;
%for i = 1:conv_inc
%    A(i) = A_0 + (A_max-A_0) * (t/t_Amax) * exp((1-t)/t_Amax);
%    
%    
%    t = t + pulse_tstep;
%end


figure(1)
plot(conv_tvec,alpha)
hold on
stem(conv_tvec, xMusf)

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
