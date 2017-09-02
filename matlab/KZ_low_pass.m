function [aout, t1] = KZ_low_pass(data, f1);

%t = 1:365;
%t = [t, t+365, t+365*2];
%data = reshape(data, 1, length(data));
%data = [data, data, data];

t = 1:length(data);

t1 = t;
data1 = data;
for i = 1:f1(2)
	[d1, temp] = moving_average_2(data1, f1(1));
	t1 = t1(1)-1+temp;
	data1 = d1;
end

aout = data1;
