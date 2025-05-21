function [outputArg1,outputArg2] = untitled(inputArg1,inputArg2)
% coefficient pattern for central difference method
forward_coefs = zeros(N-2,N);
for i=2:N-1
    forward_coefs(i,i-1) = -1;
    forward_coefs(i,i) = 0;
    forward_coefs(i,i+1) = 1;
end
end

