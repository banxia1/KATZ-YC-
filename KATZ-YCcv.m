function KATZ-YCcv(gamadd,gamall,beta,k)
%KATZ-YCcv(1,1,0.01,2)
%predict trait-related gene based on KATZ-YC in the term of leave-one-out cross validation
%A: Binary relations between trait and gene, 1st column:trait, 2nd column:gene
A=textread('knowntraitgeneinteraction.txt');
% nt:the number of traits
% ng:the number of gene
% pp:the number of known diseae-gene associations
nt=max(A(:,1)); 
ng=max(A(:,2));
[pp,qq]=size(A);
%interaction: adjacency matrix for the trait-gene association network
%interaction(i,j)=1 means gene j is related to trait i
for i=1:pp
        interaction(A(i,1),A(i,2))=1;
end

save interaction interaction;
%implement leave-one-out cross validation
for cv=1:pp 
    % obtain training sample
    load interaction;
    interaction(A(cv,1),A(cv,2))=0;
   
%calculate gamad for Gaussian kernel calculation
    for i=1:nt
        sd(i)=norm(interaction(i,:))^2;
    end
    gamad=nt/sum(sd')*gamadd;
    
 %calculate gamal for Gaussian kernel calculation
    for i=1:ng
        sl(i)=norm(interaction(:,i))^2;
    end
    gamal=ng/sum(sl')*gamall;
    
    %calculate Gaussian kernel for the similarity between trait: kd
    for i=1:nt
        for j=1:nt
    pkd(i,j)=exp(-gamad*(norm(interaction(i,:)-interaction(j,:)))^2);
        end
    end
    
    %calculate Gaussian kernel for the similarity between gene: km
        for i=1:ng
            for j=1:ng
                km(i,j)=exp(-gamal*(norm(interaction(:,i)-interaction(:,j)))^2);
            end
        end 

        for i=1:nt
            for j=1:nt
                kd(i,j)=1/(1+exp(-15*pkd(i,j)+log(9999)));
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


% obtain the score of tested  trait-gene interaction
finalscore=F(A(cv,2),A(cv,1));
% make the score of seed  trait-gene interactions as zero
for i=1:nt
    for j=1:ng
        if interaction(i,j)==1
           F(j,i)=-10000;
        end
    end
end

% obtain the position of tested trait-gene interaction as variable globalposition(1,cv),
[ll1,mm1]=size(find(F>=finalscore));
[ll2,mm2]=size(find(F>finalscore));
globalposition(1,cv)=ll2+1+(ll1-ll2-1)/2;

end
 save('globalposition.mat','globalposition');   
end


        
        
        
    
   



