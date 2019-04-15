function [a_out, alpha_out, difference_out] = GradientDescentLoop(y,H,X,lambda,a,alpha_initial,tolerance)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
m=0; x=X;
iterations=0; iterations_out=0;
a_out=zeros(size(a)); alpha_out=alpha_initial;
error=0; mse_old=0; mse_new=0; difference=0; count=0;
    m=X*(y-H);
m=m./(length(y));
slope_old=m;
theta_old =a;
theta_new=zeros(size(a));
alpha=alpha_initial;
a_copy=a;
if(lambda==0)
    flag=0;
else
    flag=1;
end
mse_old=MSE(y,H,flag,lambda,a);


while 1
    theta_new=theta_old-(alpha.*slope_old');
    
    a_copy=theta_new;
    H=HypothesisI(X',a_copy');
    m=0;
    
    m=X*(y-H);
m=m./(length(y));
slope_new=m;

mse_new=MSE(y,H,flag,lambda,a_copy);


   difference=abs(mse_old-mse_new);
   if(difference<1e-3)
        a_out=theta_new;
        alpha_out=alpha;
        difference_out=difference;
       break;
   end
   
if(abs(mse_old)<abs(mse_new))
    alpha=0.5*alpha;
else
    
slope_old=slope_new;
theta_old=theta_new;
mse_old=mse_new;  
end
       
       
       

       
     
error(iterations+1)=mse_old;
iterations=iterations+1;

end
% 
% if ~(iterations==0)
% figure()
% plot([1:1:iterations], error);
% end

% if ~(iterations==0)
% figure()
% plot([1:1:iterations], thet);
% end
end
