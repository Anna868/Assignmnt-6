function [HI] =HypothesisI(X,a)
%The first function in the set of four proposed Hypotheses
%
%
%It is in the standard form of linear regression with multiple variables
%
%
%i.e.: the output is of the form: H(x0,x1,...,x19)=theta0.x0+theta1.x1+...+theta19.x19
%
%
%Inputs are: 1) Matrix of scaled features 'X' 2)Vector of coefficients 'a'
HI=zeros(size(X(1,:)));
HI=X*a;

end

