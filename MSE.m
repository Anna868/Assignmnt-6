function [mse] = MSE(y,H,flag,lambda,a)
% This function computes the mean squared error of a hypothesis function H
% using the known values of the output y
mse=0; m=0;reg=0;
if~(flag==1)
for i=1:1:length(y)
    mse=mse+(y(i)-H(i))^2;
    m=m+1;
end
mse=mse./(2*m);
else
    for i=1:1:length(y)
    mse=mse+(y(i)-H(i))^2;
    m=m+1;
    end
    for i=2:1:length(a)
    reg=reg+abs(a(i));
    end
    mse=(mse./(2*m))+(lambda)*reg;
end 
        
end