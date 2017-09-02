
dT = 30; % min 
N = 7*24*60/dT; t = (0:N-1)*dT; % data for 7 days
pnoise = 0.30;
T1 = 12.4*60; T2 = 24*60; T3 = 15*24*60; Tc = 10*60; % min
xn = 5 + 3*cos(2*pi*t/T1) + 2*cos(2*pi*t/T2) + 1*cos(2*pi*t/T3);
xn = xn + pnoise*max(xn-mean(xn))*(0.5 - rand(size(xn)));   
[xs,c,h,Cx,f] = lanczosfilter(xn,dT,1/Tc,[],'low');  
subplot(211), plot(t,xn,t,xs), legend('noisy','smooth'), axis tight
subplot(212), plot(f,h,f,abs(Cx)/max(abs(Cx)),...
   [1 1]/Tc,[min(h) max(h)],'-.',...
   [1/T1 1/T2 1/T3],([1/T1 1/T2 1/T3]<=1/Tc),'o'), axis tight  
