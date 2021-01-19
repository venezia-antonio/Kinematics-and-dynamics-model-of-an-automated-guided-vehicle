
close all
clear all
clc
syms Rot theta
syms Rt
syms a14 a23
syms phi
syms M Iz
syms bfront brear bleft bright
syms i j k
syms e1 e2
syms eq1 eq2 eq3 eq4
syms t
syms Cm1 Cm2 Cm3 Cm4
syms F1 F2 F3 F4
syms Cmx
syms Cm1x Cm1y Cm1z Cm2x Cm2y Cm2z Cm3x Cm3y Cm3z Cm4x Cm4y Cm4z
syms accTheta accInerz ax_Inerz ay_Inerz
%%
%Assegnazione dati derivanti dallo studio geometrico del veicolo
Rt = 0.08;
a14 = pi/4;
a23 = 3*pi/4;
M = 600;
Iz = M/12*(0.5^2 + 1^2);
Bfront = 0.5;
Brear  = 0.5;
Bleft  = 0.25;
Bright = 0.25;
%%
%Definizione della base ortonormale
i=[1 0 0]';
j=[0 1 0]';
k=[0 0 1]';

%%
%Definizione matrice di rotazione
Rot(theta)=[cos(theta) -sin(theta)  0;
            sin(theta)  cos(theta)  0;
                0           0       1];
%%
%Definizione sistema di riferimento solidale al rullino         
e1=Rot(a14)*i;
e2=Rot(a23)*i;

%%
%Definizione vettori coppia motrice
Cm1 = [Cm1x, Cm1y, Cm1z];
Cm2 = [Cm2x, Cm2y, Cm2z];
Cm3 = [Cm3x, Cm3y, Cm3z];
Cm4 = [Cm4x, Cm4y, Cm4z];

%%
%Equazione di equilibrio alla rotazione per ogni singola ruota
eq1 = Cm1 + cross(k*Rt, F1*e1);
eq2 = Cm2 + cross(k*Rt, F2*e2);
eq3 = Cm3 + cross(k*Rt, F3*e2);
eq4 = Cm4 + cross(k*Rt, F4*e1);

%%
%Calcolo delle forze scambiate con il terreno Fi
solF=solve([eq1(1)==0,eq2(1)==0,eq3(1)==0,eq4(1)==0],[F1,F2,F3,F4]);

%%
%Definizione direzioni delle Fi
F1=solF.F1*e1;
F2=solF.F2*e2;
F3=solF.F3*e2;
F4=solF.F4*e1;
%%
%Calcolo accelerazione del centro di massa
%nel sistema di riferimento locale
 
 acc=(F1+F2+F3+F4)/M;


%%
%Calcolo accelerazione angolare del centro di massa nel sistema di
%riferimento locale
accTheta=1/Iz*(cross([-bleft;bfront;0],F1)+cross([bright;bfront;0],F2)+cross([-bleft;-brear;0],F3)+cross([bright;-brear;0],F4));
 
                   
%%
%Determinazione accelerazione,accelerazione angolare nel
%sistema di riferimento fisso

accInerz = Rot(phi)*acc;
ax_Inerz(t) = accInerz(1);
ay_Inerz(t) = accInerz(2);
accTheta_Inerz(t) = accTheta(3);

%%
%Risoluzione equazioni differenziali

%Legge oraria di rotazione Theta
syms Theta(t)
eqn=diff(Theta,t,2)== accTheta_Inerz;
DTheta = diff(Theta,t);
cond=[Theta(0)==0,DTheta(0)==0];
Theta(t)=dsolve(eqn,cond);

%Legge oraria lungo X
syms ax(t)
eqn=diff(ax,t,2)== ax_Inerz;
Dax = diff(ax,t);
cond=[ax(0)==0,Dax(0)==0];
X(t)=dsolve(eqn,cond);

%Legge oraria lungo X
syms ay(t)
eqn=diff(ay,t,2)== ay_Inerz;
Day = diff(ay,t);
cond=[ay(0)==0,Day(0)==0];
Y(t)=dsolve(eqn,cond);

%%
%Definizione vettore delle coppie
coppie=[Cm1x Cm2x Cm3x Cm4x];

%Scelta tra traiettorie predefinite
scelta = input('Inserisci:\n 1 - Traiettoria rettilinea\n 2 - Traiettoria obliqua\n 3 - Rotazione intorno al baricentro\n 4 - Traiettoria casuale -->'); 
switch scelta
    case 1
        val_coppie=[10 10 10 10];
    case 2
        val_coppie=[10 0 0 10];
    case 3
        val_coppie=[-10 10 -10 10];
    case 4
        val_coppie=[sin(t) -cos(t) -cos(t) sin(t)];
end

%Calcolo di Theta, necessario per la sostituzione simbolica
%per le leggi orarie lungo X e Y

var_ = [Rt,a14,a23,M,Iz,bfront,brear,bleft,bright];
val_ = [0.08 pi/4 3*pi/4 600 53.125 0.5 0.5 0.25 0.25];
Theta_sub(t)=subs((subs(Theta(t),coppie,val_coppie)),var_,val_);

%Sostituzione dei valori per le rispettive variabili simboliche
var = [Rt,a14,a23,M,Iz,bfront,brear,bleft,bright,phi];
val = [0.08 pi/4 3*pi/4 600 53.125 0.5 0.5 0.25 0.25,Theta_sub(t)];

X_sub(t)=subs((subs(X(t),coppie,val_coppie)),var,val);
Y_sub(t)=subs((subs(Y(t),coppie,val_coppie)),var,val);
ax_sub(t)=subs(subs(ax_Inerz(t),coppie,val_coppie),var,val);
ay_sub(t)=subs(subs(ay_Inerz(t),coppie,val_coppie),var,val);
accTheta_sub(t)=subs(subs(accTheta_Inerz(t),coppie,val_coppie),var,val);

%%
%Plotting

%Plot X(t)
figure (1)
subplot(3,1,1)
set(gca,'fontsize',20)
hold on
box on
fplot(X_sub(t),[0 100],'linewidth',2.0);
xlabel('t[s]');
ylabel('X(G) [m]');

%Plot Y(t)
subplot(3,1,2)
set(gca,'fontsize',20)
hold on
box on
fplot(Y_sub(t),[0 100],'linewidth',2.0);
xlabel('t[s]');
ylabel('Y(G) [m]');

%Plot Theta(t)
subplot(3,1,3)
set(gca,'fontsize',20)
hold on
box on
fplot(Theta_sub(t),[0 100],'linewidth',2.0);
xlabel('t[s]');
ylabel('Theta [rad]');

%Plot Traiettoria
figure (4)
set(gca,'fontsize',20)
hold on
box on
fplot(X_sub(t),Y_sub(t),[0 100],'linewidth',2.0);
xlabel('X(G) [m]');
ylabel('Y(G) [m]');

%Plot ax(G)
figure (5)
subplot(3,1,1)
set(gca,'fontsize',20)
hold on
box on
fplot(ax_sub(t),[0 100],'linewidth',2.0);
xlabel('t[s]');
ylabel('ax(t) [m/s^2]');

%Plot ay(G)
subplot(3,1,2)
set(gca,'fontsize',20)
hold on
box on
fplot(ay_sub(t),[0 100],'linewidth',2.0);
xlabel('t[s]');
ylabel('ay(t) [m/s^2]');

%Plot accTheta(G)
subplot(3,1,3)
set(gca,'fontsize',20)
hold on
box on
fplot(accTheta_sub(t),[0 100],'linewidth',2.0);
xlabel('t[s]');
ylabel('accTheta(t) [m/s^2]');


