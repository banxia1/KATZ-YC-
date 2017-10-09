function overallauc=positiontooverallauc5cv()
overallauc=zeros(3,1)
for k2=2:1:4
   KATZ-YC5cv(1,1,0.01,k2)
load position.mat;
load interaction;
[n,m]=size(interaction);
sID=textread('knowntraitgeneinteraction.txt');
[pp,qq]=size(sID);

for i=1:pp
    if i<floor(pp/5)*4+1

knowndt(i)=pp-floor(pp/5);
        unknown(i)=n*m-knowndt(i);
    else knowndt(i)=floor(pp/5)*4;
         unknown(i)=n*m-knowndt(i);
    end
end

for i=1:pp
    for j=1:m*n
        if j<unknown(i)+1
            rank(i,j)=j;
        else rank(i,j)=unknown(i)+1;
        end
    end
end

for k=1:m*n-floor(pp/5)*4
    tp=0;
    for t=1:pp
        if position(1,t)<=k
            tp=tp+1;
        end
    end
    tpr(1,k)=tp/pp;
    if k<m*n-pp+floor(pp/5)+1
    fp=k*pp-tp;
    else fp=floor(pp/5)*4*(m*n-pp+floor(pp/5))+(pp-floor(pp/5)*4)*k-tp;
    end

     fpr(1,k)=fp/(floor(pp/5)*4*(m*n-pp+floor(pp/5)-1)+(pp-floor(pp/5)*4)*(m*n-floor(pp/5)*4-1));
end
  

    
plot(fpr,tpr),xlabel('FPR'),ylabel('TPR')
hold on
clear area;
area(1,1)=tpr(1,1)*fpr(1,1)/2;
for k=2:m*n-floor(pp/5)*4
    area(1,k)=[tpr(1,k-1)+tpr(1,k)]*[fpr(1,k)-fpr(1,k-1)]/2;
end
overallauc(k2-1,1)=sum(area);
end
legend('k=2','k=3','k=4')
end
          


