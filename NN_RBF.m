% *************************************************************************
% -------------------------------------------------------------------------
%
%   File        NN_RBF.m
%   Purpose     Radial Basis Function Network
%   Author      Vynokurova O.
%   Modified    27/08/2017 - release
%
% -------------------------------------------------------------------------
% *************************************************************************

clear all; close all; clc;

% ��������� �������� ���� ��� ����� ��������� ������� ��� � �����
load data_sunspot.dat % Monthly smoothed total sunspot number http://www.sidc.be/silso/datafiles#total
x=data_sunspot;

% ����������� �����
%x=(x-min(x))/(max(x)-min(x));        % �� 0 �� 1
x=2*((x-min(x))/(max(x)-min(x)))-1;   % �� -1 �� 1
% x=(x-mean(x))/std(x);               % ����������� �� ������ �� 
                                      % ����������� ���������

% ����������� ���������
Bias=1;                               % ���� (or 0)        
n=2;                                  % ����� �����

% ���������� ������ ����� � ��������� �������� �������� �� 1 ���� 
for t=(n+1):length(x) 
    for i=1:n+1
        Data(t,i)=x(t-(n+1-i)); 
    end;
end;

% ���������� ��������� �� ������� ������
TTotal=size(Data,1);                 % ��������� ������
TTrn=round(2/3*TTotal);              % ��������� ��������� ������
TChk=TTotal-TTrn;                    % ��������� ������� ������

trnData=Data(1:TTrn,:);              % ��������� ������
chkData=Data(TTrn+1:TChk+TTrn,:);    % ������� ������

%��������� ���������� ������� ������ ��������-�������� �������
%������ 1 - ��������� ������ ��������� ��������
maxx=max(Data(:,1));                   % ���������� ������������ ��������
minx=min(Data(:,1));                   % ���������� ��������� ��������

% ���������� ������� ������ ��������� �������� ��� n=2
d=4; %������� d �� d 
cs=linspace(minx,maxx,d);
C=[];
for i=1:d
    for j=1:d
        C=[C; cs(i) cs(j)];
    end
end
h=d^2;

%������ 2 - ��������� ������ �� ��������� ������� �������������
% [C,S] = subclust(trnData,0.4);
% h=size(C,1);                           % ������� ��������-�������� �������
% C=C(:,1:n);                            % ������� ������  

%��������� ���������� ������� ����� ��������-�������� �������
sigma_init(1:n)=4;
for i=1:h
     Sigma(i)={diag(sigma_init)}; % ��������� � ��������� ����� 
end; 

% ����������� ���������� ��������� ��������� �������� ����������� �������
% �����������
W=zeros(h+1,1); 
P=10000*eye(h+1);
alpha_w=0.99;
AddW=zeros(h+1,1);

% ����������� ���������� ��������� ��������� �������� ������ �� ����� ���
rC(1:h)=100; 
rS(1:h)=100; 
lambda_c=0.99;
lambda_s=0.99;
alpha_c=0.99;
alpha_s=0.99;
AddC=zeros(h,n);
AddSigma(1:h)={zeros(n,n)};

%���������� �����c� ��������
num_eppoch=10;                                            % ʳ������ ���� ��������
for epoch=1:num_eppoch
    for k=n+1:TTrn
        Xin=trnData(k,1:n);                               % ������� ������ �����
        [y(k), Phi, Phi_C, Phi_S]=netout(Xin',C,Sigma,W); % ���������� ������ �����
        d(k)=trnData(k,n+1);                              % ���������� ������
        e(k)=d(k)-y(k);                                   % ���������� ������� ��������
        %���������� ��������� ��������� ��������
        % ��� ������� �����������
        P = P - (P * Phi * Phi' * P) / (1 + Phi' * P * Phi);
        AddW=(P*e(k)*Phi)/(1+Phi'*P*Phi);
        % ��� ������ ��������-�������� �������  
        for i=1:h
            xc=Xin'-C(i,:)';
            J_c=2*W(i)*Phi_C(i)*Sigma{i}*xc;
            AddC(i,:)=-lambda_c*e(k)*J_c/rC(i);
            rC(i)=alpha_c*rC(i)+(J_c'*J_c);
        end;
        % ��� ����� ��������-�������� �������    
        for i=1:h
            xc=Xin'-C(i,:)';
            J_s=W(i)*Phi_S(i)*xc*xc';
            AddSigma(i)={(lambda_s*e(k)*J_s)/rS(i)};
            rS(i)=alpha_s*rS(i)+trace(J_s'*J_s);
        end;
        % ��������� ������� �����������
        W=W+AddW;
        % ��������� ������ ��������-������� �������
        C=C-AddC;
        % ��������� ����� ��������-������� �������
        for i=1:h
            Sigma(i)={Sigma{i}-AddSigma{i}};
        end;   
   end;
       % ���������� ������������������� ������� ��������
       RMSEtrn(epoch)=norm(e(1:TTrn))/sqrt(TTrn);
       % ��������� ������������������� ������� �������� � �������
       disp(['Epoch: ', num2str(epoch), ', RMSEtrn=', num2str(RMSEtrn(epoch))]);      
end;

%���������� ��������� ����������
for k=1:TChk
    Xin=chkData(k,1:n);                 % ���������� �������� �������
    y_chk(k)=netout(Xin',C,Sigma,W);    % ���������� ������ �����  
    d=chkData(k,n+1);                   % ���������� ������ 
    e_chk(k)=d-y_chk(k);                % ���������� ������� ����������
end;

% ³��������� ������� �������������
figure;
plot(1:TTotal,Data(:,n+1),'--');
hold on;
plot(1:TTrn,y(1:TTrn));
hold on;
plot(TTrn+1:TTotal,y_chk(1:TChk),'g');
title(['RMSEtrn=', num2str(RMSEtrn(epoch)), '   RMSEchk=' num2str(norm(e_chk(1:TChk))/sqrt(TChk))]);
axis tight;
figure
plot(1:TChk,chkData(1:TChk,n+1),1:TChk,y_chk(1:TChk));
title(['RMSEchk=' num2str(norm(e_chk(1:TChk))/sqrt(TChk))]);
axis tight;

% ���������� ������� ������������� �� ������� ������ 
RMSE=(norm(e_chk(1:TChk))/sqrt(TChk));                                  % Root-mean-square error
SMAPE=100/TChk*(sum((abs(e_chk))./((abs(chkData(:,n+1)')+abs(y_chk))))) % Symmetric mean absolute percentage error
