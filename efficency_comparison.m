m=50;
q1s=rand(1,1000)*0.1;
qns=rand(1,1000)*2+0.1;
tuns=rand(1,1000)*2-1;
build=zeros(1,1000);
my=zeros(1,1000);
for k=1:length(q1)
n=m;
q1=q1s(k);qn=qns(k);tun=tuns(k);lambda=qn*1.5;
nn=1:n;
if abs(tun)<10^-3
    qq_real=q1.^((n-nn)/(n-1)).*qn.^((nn-1)/(n-1));
else
    qq_real=(q1^tun*(n-nn)/(n-1)+qn^tun*(nn-1)/(n-1)).^(1/tun);
end
yi=coxian(5000,lambda,qq_real);
yi=sort(yi);
y=[yi,ones(5000,1)];
tic
objfun1_build(m,q1,qn,lam,tun,y);
build(k)=toc;
tic
objfun1(m,q1,qn,lam,tun,y,10^-10);
my(k)=toc;
end
h = figure;
histogram(build)
hold on
histogram(my)
legend({'expm','our algorithm'},'Location','northeast')
xlabel('seconds')
title('required time to evaluate likelihood')
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'efficiency_compare','-dpdf','-r0')
