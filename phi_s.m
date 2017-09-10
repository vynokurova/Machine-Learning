function y = phi_s(x,c,sigma)

y=(norm_is(x-c,sigma))*exp(-(norm_is(x-c,sigma))^2/2);