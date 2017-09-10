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

% Генерація часового ряду або можна загрузити часовий ряд з файлу
load data_sunspot.dat % Monthly smoothed total sunspot number http://www.sidc.be/silso/datafiles#total
x=data_sunspot;

% Нормалізація даних
%x=(x-min(x))/(max(x)-min(x));        % від 0 до 1
x=2*((x-min(x))/(max(x)-min(x)))-1;   % від -1 до 1
% x=(x-mean(x))/std(x);               % нормалізація на середнє та 
                                      % квадратичне відхилення

% Ініціалізація параметрів
Bias=1;                               % зсув (or 0)        
n=2;                                  % число входів

% Формування вибірки даних з допомогою елементів затримки на 1 такт 
for t=(n+1):length(x) 
    for i=1:n+1
        Data(t,i)=x(t-(n+1-i)); 
    end;
end;

% Формування навчальної та тестової вибірки
TTotal=size(Data,1);                 % Розмірність вибірки
TTrn=round(2/3*TTotal);              % Розмірність навчальної вибірки
TChk=TTotal-TTrn;                    % Розмірність тестової вибірки

trnData=Data(1:TTrn,:);              % Навчальна вибірка
chkData=Data(TTrn+1:TChk+TTrn,:);    % Тестова вибірка

%Генерація початкових значень центрів радіально-базисних функцій
%Варіант 1 - Розміщення центрів рівномірною решіткою
maxx=max(Data(:,1));                   % обчислення максимальних елементів
minx=min(Data(:,1));                   % обчислення мінімальних елементів

% Формування матриці центрів рівномірною решіткою для n=2
d=4; %решітка d на d 
cs=linspace(minx,maxx,d);
C=[];
for i=1:d
    for j=1:d
        C=[C; cs(i) cs(j)];
    end
end
h=d^2;

%Варіант 2 - Розміщення центрів за допомогою функції кластеризації
% [C,S] = subclust(trnData,0.4);
% h=size(C,1);                           % кількість радіально-базисних функцій
% C=C(:,1:n);                            % матриця центрів  

%Генерація початкових значень ширин радіально-базисних функцій
sigma_init(1:n)=4;
for i=1:h
     Sigma(i)={diag(sigma_init)}; % структура з матрицями ширин 
end; 

% Інціалізація початкових параметрів алгоритма навчання синаптичних вагових
% коефіцієнтів
W=zeros(h+1,1); 
P=10000*eye(h+1);
alpha_w=0.99;
AddW=zeros(h+1,1);

% Інціалізація початкових параметрів алгоритмів навчання центрів та ширин РБФ
rC(1:h)=100; 
rS(1:h)=100; 
lambda_c=0.99;
lambda_s=0.99;
alpha_c=0.99;
alpha_s=0.99;
AddC=zeros(h,n);
AddSigma(1:h)={zeros(n,n)};

%Формування процеcу навчання
num_eppoch=10;                                            % Кількість епох навчання
for epoch=1:num_eppoch
    for k=n+1:TTrn
        Xin=trnData(k,1:n);                               % Вхідний вектор даних
        [y(k), Phi, Phi_C, Phi_S]=netout(Xin',C,Sigma,W); % Обчислення виходу мережі
        d(k)=trnData(k,n+1);                              % Навчальний вектор
        e(k)=d(k)-y(k);                                   % Обчислення похибки навчання
        %Розрахунок параметрів алгоритмів навчання
        % Для вагових коефіцієнтів
        P = P - (P * Phi * Phi' * P) / (1 + Phi' * P * Phi);
        AddW=(P*e(k)*Phi)/(1+Phi'*P*Phi);
        % Для центрів радіально-базисних функцій  
        for i=1:h
            xc=Xin'-C(i,:)';
            J_c=2*W(i)*Phi_C(i)*Sigma{i}*xc;
            AddC(i,:)=-lambda_c*e(k)*J_c/rC(i);
            rC(i)=alpha_c*rC(i)+(J_c'*J_c);
        end;
        % Для ширин радіально-базисних функцій    
        for i=1:h
            xc=Xin'-C(i,:)';
            J_s=W(i)*Phi_S(i)*xc*xc';
            AddSigma(i)={(lambda_s*e(k)*J_s)/rS(i)};
            rS(i)=alpha_s*rS(i)+trace(J_s'*J_s);
        end;
        % Оновлення вагових коефіцієнтів
        W=W+AddW;
        % Оновлення центрів радіально-базисної функції
        C=C-AddC;
        % Оновлення ширин радіально-базисної функції
        for i=1:h
            Sigma(i)={Sigma{i}-AddSigma{i}};
        end;   
   end;
       % Обчислення середньоквадратичної похибки навчання
       RMSEtrn(epoch)=norm(e(1:TTrn))/sqrt(TTrn);
       % Виведення середньоквадратичної похибки навчання в консоль
       disp(['Epoch: ', num2str(epoch), ', RMSEtrn=', num2str(RMSEtrn(epoch))]);      
end;

%Формування процедури тестування
for k=1:TChk
    Xin=chkData(k,1:n);                 % формування вхідного вектору
    y_chk(k)=netout(Xin',C,Sigma,W);    % Обчислення виходу мережі  
    d=chkData(k,n+1);                   % Навчальний сигнал 
    e_chk(k)=d-y_chk(k);                % Обчислення похибки тестування
end;

% Візуалізація процесу прогнозування
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

% Обчислення критеріїв прогнозування на тестовій вибірці 
RMSE=(norm(e_chk(1:TChk))/sqrt(TChk));                                  % Root-mean-square error
SMAPE=100/TChk*(sum((abs(e_chk))./((abs(chkData(:,n+1)')+abs(y_chk))))) % Symmetric mean absolute percentage error
