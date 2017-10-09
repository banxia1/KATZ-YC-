function KATZ-YC2cv(gamadd,gamall,beta,k)
%KATZ-YC2cv(1,1,0.01,2)
%predict trait-related gene based on KATZ-YC  in the term of 2-fold cross validation
%A: Binary relations between trait and gene, 1st column:trait, 2nd column:gene
A=textread('knowntraitgeneinteraction.txt');
% nd:the number of traits
% nm:the number of gene
% pp:the number of known trait-gene associations
nd=max(A(:,1)); 
nm=max(A(:,2));
[pp,qq]=size(A);
%interaction: adjacency matrix for the trait-gene association network
%interaction(i,j)=1 means gene j is related to trait i
for i=1:pp
        interaction(A(i,1),A(i,2))=1;
end

save interaction interaction;



%implement 2-fold cross validation
x=randperm(pp)';
T=1;

for cv=1:2
load interaction interaction;
    if cv<2
        B=A(x((cv-1)*floor(pp/2)+1:floor(pp/2)*cv),:);
% obtain training sample
for i=1:floor(pp/2)
        interaction(B(i,1),B(i,2))=0;
    end
    else B=A(x((cv-1)*floor(pp/2)+1:pp),:);
        % obtain training sample
for i=1:pp-floor(pp/2)
        interaction(B(i,1),B(i,2))=0;
    end
    end
    
   
%calculate gamad for Gaussian kernel calculation
    for i=1:nd
        sd(i)=norm(interaction(i,:))^2;
    end
    gamad=nd/sum(sd')*gamadd;
    
 %calculate gamal for Gaussian kernel calculation
        for i=1:nm
        sl(i)=norm(interaction(:,i))^2;
    end
    gamal=nm/sum(sl')*gamall;
    
    %calculate Gaussian kernel for the similarity between trait: kd
    for i=1:nd
        for j=1:nd
    pkd(i,j)=exp(-gamad*(norm(interaction(i,:)-interaction(j,:)))^2);
        end
    end
    
    %calculate Gaussian kernel for the similarity between gene: km
        for i=1:nm
        for j=1:nm
    km(i,j)=exp(-gamal*(norm(interaction(:,i)-interaction(:,j)))^2);
        end
        end 

        for i=1:nd
for j=1:nd
    kd(i,j)=1/(1+exp(-12*pkd(i,j)+log(9999)));
end
        end

if k==2
F=beta*interaction'+(beta^2)*(km*interaction'+interaction'*kd);
else if k==3
F=beta*interaction'+(beta^2)*(km*interaction'+interaction'*kd)+(beta^3)*(interaction'*interaction*interaction'+km*km*interaction'+km*interaction'*kd+interaction'*kd*kd);
else if k==4
F=beta*interaction'+(beta^2)*(km*interaction'+interaction'*kd)+(beta^3)*(interaction'*interaction*interaction'+km*km*interaction'+km*interaction'*kd+interaction'*kd*kd)+(beta^4)*(km*km*km*interaction'+interaction'*interaction*km*interaction'+km*interaction'*interaction*interaction'+interaction'*kd*interaction*interaction')+(beta^4)*(interaction'*interaction*interaction'*kd+km*km*interaction'*kd+km*interaction'*kd*kd+interaction'*kd*kd*kd);
end
end
end

F=F';

[size1B,size2B]=size(B);
% obtain the score of tested  trait-gene interaction
for i=1:size1B
finalscore(i,1)=F(B(i,1),B(i,2));
end
% make the score of seed  trait-gene interactions as zero
for i=1:nd
    for j=1:nm
        if interaction(i,j)==1
           F(i,j)=-10000;
        end
    end
end


for qq=1:size1B
% obtain the position of tested trait-gene interaction as variable position(1,cv), 
[ll1,mm1]=size(find(F>=finalscore(qq)));
[ll2,mm2]=size(find(F>finalscore(qq)));
position(1,T)=ll2+1+(ll1-ll2-1)/2;
T=T+1;
end

end
save('position.mat','position');  
end



        
        
        
    
   



