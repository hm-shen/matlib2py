function [moveav, tnew] = moving_average_2(d, t)

% this is only for simple moving average

len = length(d);
beginno = ceil(t/2);
endno = len - floor(t/2);

moveav = 0;

tnew = beginno : endno ;

for i = 1:t
	moveav = moveav + d(1+i-1:len-t+i);
end

	moveav = moveav/t;
