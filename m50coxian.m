rng(795)
b=0.015;
c=0.01;
m=50;
nn=1:m;
qq_real=b+c*nn;
lambda=1.6; 
% 100 bootstrap for estimator
boot=zeros(100,9);
boot_time=zeros(size(boot,1),1);
parfor bt=1:size(boot,1)
    %simulate data
yi=coxian(1000,lambda,qq_real);
% yi=sort(yi);
% y=[yi,ones(5000,2)]; 
[a,b]=hist(yi, unique(round(yi,2))); y=[b,a']; 
% tic
dx=[];
for mm=3:7
      dx=[dx;estimatefun_fix(10*mm,y),10*mm];
end
[~,i]=min(dx(:,1));
% boot_time(bt)=toc;
%record the each bootstrap result
 boot(bt,:)=dx(i,:);   
end

