function y = phi_c(x,c,sigma)

y=(norm_is(x-c,sigma))*exp(-(norm_is(x-c,sigma))^2/2);