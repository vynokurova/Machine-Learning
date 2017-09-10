function [y, Phi, Phi_C, Phi_S]=netout(x,C,Sigma,W)
% Обчислення виходу мережі та компонентів алгоритмів навчання

[h,n]=size(C);
for i=1:h
    Phi(i,1)=phi(x,C(i,:)',Sigma{i});
end;
Phi(h+1,1)=1;

y=W'*Phi;

for i=1:h
    Phi_C(i,1)=phi_c(x,C(i,:)',Sigma{i});
end;

for i=1:h
    Phi_S(i,1)=phi_s(x,C(i,:)',Sigma{i});
end;