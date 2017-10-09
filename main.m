for cv=1:100
KATZ-YC5cv(1,1,0.01);
overallauc(cv)=positiontooverallauc();
end
save overallauc overallauc
a=mean(overallauc);
b=std(overallauc);