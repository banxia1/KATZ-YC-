function overallauc=positiontooverallauc()
overallauc=zeros(3,1)
for k1 = 2:1:4
KATZ-YCcv(1,1,0.01,k1)
load globalposition.mat;
load interaction;
[n,m]=size(interaction);
sID=textread('knowntraitgeneinteraction.txt');
[pp,qq]=size(sID);


for i=1:pp
if globalposition(i)>m*n-pp+1
globalposition(i)=m*n-pp+1;
end
end
for k=1:m*n-pp+1
    tp=0;
    for t=1:pp
        if globalposition(1,t)<=k
            tp=tp+1;
        end
    end
    tpr(1,k)=tp/pp;
    fp=k*pp-tp;
    fpr(1,k)=fp/(pp*(m*n-pp));
end
plot(fpr,tpr),xlabel('FPR'),ylabel('TPR')
hold on
clear area;
area(1,1)=tpr(1,1)*fpr(1,1)/2;
for k=2:m*n-pp+1
    area(1,k)=[tpr(1,k-1)+tpr(1,k)]*[fpr(1,k)-fpr(1,k-1)]/2;
end
overallauc(k1-1,1)=sum(area);
end
legend('k=2','k=3','k=4')
end

          
