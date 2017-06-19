% test convolution
clear

MF_firing = [0.
0.
0.
0.
1.
0.
0.
1.
0.
0.
0.
0.
0.
0.
0.
1.
0.
0.
0.
0.
0.
0.
0.
0.
0.
0.
0.
0.
0.
1.
0.
0.
1.
0.
0.
0.
0.
0.
0.
1.
0.
0.
1.
0.
0.
0.
0.
0.
0.
0.
1.
0.
0.
0.
0.
0.
0.
0.
0.
0.
0.
0.
0.
0.
1.
0.
0.
1.
0.
0.
1.
0.
0.
1.
0.];

% fpulse
t = 0;
Pprime = 5.0;
Tcprime = 4E-9;
t=0;
for j = 1:125
      f_pulse(j) = t/Tcprime * exp(1-t/Tcprime); 
      t = t+1E-9;
end 

%for i = 1:length(a)
%    MF_firing(i) = a(i);
%end

%MF_firing(1) = 1;
%MF_firing(10)= 1;

fconv = conv(MF_firing,f_pulse);

u = f_pulse;
v = MF_firing;

fconv2 = conv(v,u);

k = 1;
w(length(u)+length(v)-1)=zeros;
% pad zeros
v(length(v)+1:length(w))=zeros;
for i = 1:size(w,2)
    w(k) = 0;
    for j = 1:size(u,2)
        if (k-j+1) >= 1
            w(k) = w(k) + u(j)*v(k-j+1);
        end
    end
    k = k + 1;
end



plot(w/max(w))
%hold on
%plot(f_pulse)