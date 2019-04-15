function f = Gauss(mu,sigma,x)
f=((2*pi)^(1/2))*sigma;
f=f.^-1;
mat=-(x-mu).^2;
mat=mat./(2*sigma^2);
f=f.*(exp(mat));
end

