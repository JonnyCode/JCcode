function [ReturnVal] = factorial(InputVal)

%	[ReturnVal] = factorial(InputVal)


ReturnVal = 1;

for n=1:InputVal
	ReturnVal = ReturnVal * n;
end
