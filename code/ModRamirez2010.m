% Define variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fl: value of the force-length relationship at a given stretch
% fV: value of the force-length relationship at a given stretch
% concA: concentration of species A
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% main function
function [f_total] = ModRamirez2010(ref_signal)
close all;

%% length of muscle at the instant of the isometric contraction
lam=1.4;
[fl] = forceLength(lam);

%% voltage
V = 4.;
% disregard the "external activation"
fV = 1.;

%% time array
t = 0:0.001:1.2;
%% single puls
[fpulse,time] = forceFrequency(t);
f_total = fl*fV*fpulse;

figure,
plot(time,f_total)
xlabel('time [s]'), ylabel('force = f_\lambda \Sigma_i \delta_i f_p_u_l_s_e')
hold on
%plot(time,ref_signal)

f_total_max = max(f_total);
f_total_norm = f_total/f_total_max;
figure,
plot(time,f_total_norm)
xlabel('time [s]'), ylabel('normalized force = (f_\lambda \Sigma_i \delta_i f_p_u_l_s_e) / f_max')
end

%% force-LENGTH relationship
function [fl] = forceLength(lam)
[lamOpt,beta] = const_fl();

fl = exp(-(lam-lamOpt).^2/(2.*(1-beta)^2));
end

function [f_total, t_prime] = forceFrequency(t)
[Pprime,Tc,fr,r,c] = const_ft();
% fpulse = Pprime*t/Tc*exp(1.-t/Tc);

%% Motor unit recruitment

% secind batch of firings w/ randomized n's
% n_MU1 = randi([0 1],size(t,1),size(t,2));
% n_MU2 = randi([0 1],size(t,1),size(t,2));
% n_MU3 = randi([0 1],size(t,1),size(t,2));

% third batch of firings w/ firings from Leonardo's file ...mat
load 255MiLosa150414150101rampa20p1_decomp.mat

% maxsize_MUfiring = cellfun('size',MUPulses);
size_MU1 = size(MUPulses{1}); size_MU2 = size(MUPulses{2}); size_MU3 = size(MUPulses{3}); size_MU4 = size(MUPulses{4}); size_MU5 = size(MUPulses{5});
size_MU = [size_MU1,size_MU2,size_MU3,size_MU4,size_MU5];

for nMU = 1:size(MUPulses,2)
    eval(['n_MU', int2str(nMU),'(MUPulses{nMU}) = 1;']);
end


csvwrite('mufirings.txt', n_MU1')

%% twitch function
% wider one
f_pulse_i = t/Tc.*exp(1.-t/Tc); % normalized by Pprime
f_pulse = resample(f_pulse_i,2048,10000); %resample to match EMG recording rate

%% convolution of the firing times w/ the single twitch
list_nMU = who('n_MU*');
max_l = -1;
for n = 1:size(list_nMU,1)
    temp=eval(list_nMU{n});
    eval([ 'f_MU',int2str(n),' = conv(temp,f_pulse);']);
    temp = eval([ 'f_MU',int2str(n)]);
%     max_l = max(max_l,length(temp));
end
max_l = max(max_l,length(ref_signal)); % max length of the array such that the others are to be accodingly arranged

%% append the f_MU* arrays with zeros to make all the array lengths equal - mainly for plotting purposes
list_fMU = who('f_MU*');
for n = 1:size(list_fMU)
    eval(['f_MU',int2str(n),'_sameSize=[f_MU',int2str(n),', zeros(1,max_l-length(f_MU',int2str(n),'))];']);
end

%% sum up force contribution from every MU to get the total force
f_total = zeros(size(f_MU1_sameSize));
for n = 1:size(list_fMU)
    eval(['f_total = f_MU',int2str(n),'_sameSize + f_total ;']);
end

%% plots
t_f_total = linspace(0,1.2,size(f_total,2));
t_ref_signal= linspace(0,1.2,size(ref_signal,2));
t_prime = t_f_total;
%% increase font size on plots
set(groot, 'defaultAxesFontSize',14)
set(groot,'defaultTextFontSize',14) 

% hold on
close all;
figure
plot(t_f_total,f_total,':','LineWidth',0.7)
hold on
plot(t_ref_signal,ref_signal/max(ref_signal)*max(f_total),'LineWidth',2)
hold on
for n = 1:size(list_fMU)
    eval(['plot(t_f_total, f_MU',int2str(n),'_sameSize)']);
    hold on
end
legend('Total','Ref. signal','MU1','MU2','MU3','MU4','MU5','MU6','MU7','MU8','MU9','MU10','MU11','MU12')
xlabel('time [s]','FontSize',14), ylabel('force = \Sigma \delta_i f_p_u_l_s_e')

end

%% Constants
function [lamOpt,beta] = const_fl()
%% fl
lamOpt = 1.;
beta = .83616;
end

function [Pprime,Tc,fr,r,c] = const_ft()
%% f-t
Pprime = 0.110; % [N]
Tc = 0.4; % [s]
% Tc = 0.4; % [s]
fr = 90; % [Hz]
r = 1.0535; % [-]
c = 1.1245; % [-]
end
