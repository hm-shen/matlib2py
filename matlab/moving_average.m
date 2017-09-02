function [moveav, tnew] = moving_average(d, t)

% this is only for simple moving average
d = reshape(d, length(d), 1);
%t = reshape(t, length(t), 1);

%t = length(d);
len = length(d);
beginno = ceil(t/2);
endno = len - floor(t/2);

moveav = [];
xtemp = 0;

tnew = beginno : endno ;

for i = beginno : endno

	if mod(t, 2) == 1

	xtemp = xtemp + sum(d(i - floor(t/2) : i - floor(t/2) + t - 1));
	
	else 

        xtemp = xtemp + sum(d(i - floor(t/2) + 1 : i - floor(t/2) + t));

	end

	
	moveav = [moveav; xtemp/t];
	
	xtemp = 0; 
end

