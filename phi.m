function y = phi(x,c,sigma)
% Обчислення радіально-базисної функції активації
y=exp(-(norm_is(x-c,sigma))^2/2);