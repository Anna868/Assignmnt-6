clearvars all;
clear all;
data_set_complete=readtable('house_data_complete.csv','ReadVariableNames',false);
data_set_complete=table2array(data_set_complete(:,3:end));
 s=size(data_set_complete);

 %Extracting the 19 features as individual vectors from data set matrix and
 %scaling them via Normalization

 y_complete=(data_set_complete(:,1));
 y_complete=(y_complete-mean(y_complete))./std(y_complete);
 x0=(ones(size(y_complete)));
 s=size(data_set_complete);
 X=zeros(s(1),s(2)-1);
 for j=2:1:size(X,2)+1
     temp=(data_set_complete(:,j));
     x=(temp-mean(temp))./std(temp);
     X(:,j-1)=x';
 end
 X=X';


 %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Corr_x=corr(X');
% corrplot(int32(X'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cov_x=cov(X');
[U,S,V]=svd(Cov_x);
lambda_mat=inv(U)*Cov_x*U;
for i=1:size(Cov_x,1)
    lambda(i)=lambda_mat(i,i);
end

for i=1:length(lambda)
    alpha(i)=(sum(lambda(1:i))/sum(lambda));
end

figure()
subplot(1,2,1)
plot([1:length(lambda)],alpha, '-o'), xlabel('Eigenvalue Index'), ylabel('Variance Percent')
hold on
plot([1:length(lambda)],ones(length(lambda))*(1-0.001),'r')
subplot(1,2,2)
plot([1:length(lambda)],alpha, '-o'), xlabel('Eigenvalue Index'), ylabel('Variance Percent')
hold on
plot([1:length(lambda)],ones(length(lambda))*(1-0.01),'r')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k=16;
z=U*X;
z=z(1:k,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xapprox=U(1:k,1:k)*z;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 for i=1:k
     for j=1:size(X,2)
     err=sum((Xapprox(i,j)-z(i,j)).^2)./(k*size(z,2));
     end
 end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z0=ones(1,size(z,2));
Z=zeros(size(z,1)+1,size(z,2));
Z(1,:)=z0;
Z(2:end,:)=z;
aI=zeros(1,k+1);
HI=Z'*aI';
alpha_vector=[1:-0.1:0.1];
 for i=1:1:length(alpha_vector)
    [a_out,alpha_out, difference]=GradientDescentLoop(y_complete, HI,Z,0, aI,alpha_vector(i), 5e-4);
    d(i)=difference;
    i
end
s=size(d);
for i=1:1:s(2)
    i
pos=find(abs(d(:,i))==min(abs(d(:,i))));
if(length(pos)~=1)
best_alpha=alpha_vector(pos(1));
else
    best_alpha=alpha_vector(pos);
end
end
   [a_out1,alpha_out1, difference1]=GradientDescentLoop(y_complete, HI,Z,0, aI,best_alpha, 5e-4);

  HGrad=HypothesisI(Z',a_out1');
  Model_error=MSE(y_complete,HGrad,0,0,a_out1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xc=zeros(size(X,1)+1,size(X,2));
Xc(1,:)=x0;
Xc(2:end,:)=X;
aI=zeros(1,size(Xc,1));
HI=Xc'*aI';
alpha_vector=[1:-0.1:0.1];
 for i=1:1:length(alpha_vector)
    [a_out,alpha_out, difference]=GradientDescentLoop(y_complete, HI,Xc,0, aI,alpha_vector(i), 5e-4);
    d(i)=difference;
    i
end
s=size(d);
for i=1:1:s(2)
    i
pos=find(abs(d(:,i))==min(abs(d(:,i))));
if(length(pos)~=1)
best_alpha=alpha_vector(pos(1));
else
    best_alpha=alpha_vector(pos);
end
end
   [a_out2,alpha_out2, difference1]=GradientDescentLoop(y_complete, HI,Xc,0, aI,best_alpha, 5e-4);

  HGrad2=HypothesisI(Xc',a_out2');
  Model_error2=MSE(y_complete,HGrad2,0,0,a_out2);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Kcluster=[2:5];
 for i=1:length(Kcluster)
  [idx,~,Sum]=kmeans(X',Kcluster(i));
  Sumd(i)=sum(Sum)/Kcluster(i);
    figure()
[silh,h] = silhouette(X',idx);
h = gca;
h.Children.EdgeColor = [.8 .8 1];
xlabel 'Silhouette Value'
ylabel 'Cluster'
cluster1(i) = mean(silh);
 end
 bestCluster1=Kcluster(find(Sumd==max(Sumd)));

%   bestCluster1=Kcluster(find(cluster1==max(cluster1)));
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Kcluster=[2:5];
   for i=1:length(Kcluster)
[idx,~,Sum]=kmeans(Z',Kcluster(i));
  Sumd(i)=sum(Sum)/Kcluster(i);
%   figure()
% [silh,h] = silhouette(Z',idx);
% h = gca;
% h.Children.EdgeColor = [.8 .8 1];
% xlabel 'Silhouette Value'
% ylabel 'Cluster'
% cluster2(i) = mean(silh);
   end
   bestCluster2=Kcluster(find(Sumd==max(Sumd)));
%   bestCluster2=Kcluster(find(cluster2==max(cluster2)));
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 f=ones(k,size(y_complete,1));
  for i=1:k
     mu=mean(z(i,:));
     sigma=std(z(i,:));
     ftemp=Gauss(mu,sigma,z(i,:));
     f(i,:)=ftemp';
  end

  krand=floor(size(X,1).*rand(1,1));
  index_rand=floor(size(y_complete,1).*rand(1,1));
  
  for i=1:k
     mu=mean(z(i,:));
     sigma=std(z(i,:));
     ftemp=Gauss(mu,sigma,X(krand,index_rand));
     detection(i)=ftemp;
  end
  
  eta=0.05;
  detection=(detection>eta);
  if(detection==1)
      AnomalyFlag=0;
  else
      AnomalyFlag=1;
  end