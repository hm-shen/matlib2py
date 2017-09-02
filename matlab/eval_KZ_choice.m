function cutperiod = eval_KZ_choice(len, kzchoice, str)
% this is to use a delta function to do FFT and 
% to do KZ lowpass filtered to determine the e-folding
% freq, deemed as the cut-off for low-pass filtering
% default is low pass filter setting up

if nargin < 3
	str = 'low';
end

if size(kzchoice) ~= 2
	disp('the 2nd input must be a two-element array');
	return;
end

x = zeros(len, 1);

if mod(len, 2) == 0
x(len/2+1, 1) = 1e8;
else
x(ceil(len/2), 1) = 1e8;
end

[f1, f2] = powersp(x, length(x));

[a, t] = KZ_low_pass(x, kzchoice);

if str(1) == 'h' | str(1) == 'H'
a = x(t) - a;
end

[f3, f4] = powersp(a, length(x));

(f2(1)-f2(end))/f2(1)
[t1, t2] = min(abs(f2(1)*exp(-1) - f4));

% the unit is dt: time interval 
cutperiod = len/t2;

