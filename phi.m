function y = phi(x,c,sigma)
% ���������� ��������-������� ������� ���������
y=exp(-(norm_is(x-c,sigma))^2/2);