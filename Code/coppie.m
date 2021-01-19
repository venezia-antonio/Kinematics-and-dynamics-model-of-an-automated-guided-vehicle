close all
clear all
clc
syms Rt Rr Rm 
syms a14 a23
syms M Iz
syms bfront brear bleft bright
syms i j k
syms e1 e2
syms Rot theta
syms vGx vGy vG omega omegaz
syms omegaR1 omegaR2 omegaR3 omegaR4
syms omegaR11 omegaR22 omegaR32 omegaR41
syms omegaM1 omegaM2 omegaM3 omegaM4 omegaM1x(t) omegaM2x(t) omegaM3x(t) omegaM4x
syms eq1 eq2 eq3 eq4
syms vGloc omegaloc aGloc
syms t
syms Cm1 Cm2 Cm3 Cm4
syms F1 F2 F3 F4
syms Cm1x Cm1y Cm1z Cm2x Cm2y Cm2z Cm3x Cm3y Cm3z Cm4x Cm4y Cm4z
syms acc accang
syms axGloc ayGloc omegapuntoInerz
syms aGinerz
syms solA solB solC solD
Rm=Rt-Rr;
%%
%Definizione base ortonormale
i=[1 0 0]';
j=[0 1 0]';
k=[0 0 1]';

%%
%Definizione matrice di rotazione
Rot(theta)=[cos(theta) -sin(theta) 0;
            sin(theta)  cos(theta) 0;
            0               0      1];
%%
%Rotazione del sistema di riferimento solidale al rullino
e1=Rot(a14)*i;
e2=Rot(a23)*i;

%%
%Definizione equazioni cinematiche
vG=[vGx;vGy;0];
omega=[0;0;omegaz];
rm=[0;0;Rm];
r=[0;0;Rr];

omegaR1 = omegaR11*e1;
omegaR2 = omegaR22*e2;
omegaR3 = omegaR32*e2;
omegaR4 = omegaR41*e1;
omegaM1 = omegaM1x(t)*i; 
omegaM2 = omegaM2x(t)*i ;
omegaM3 = omegaM3x(t)*i ;
omegaM4 = omegaM4x*i; 

eq1= cross((omegaR1 + omegaM1 + omega),r) + cross((omegaM1 + omega),rm) - vG - cross(omega,[-bleft;bfront;0]);           
eq2= cross((omegaR2 + omegaM2 + omega),r) + cross((omegaM2 + omega),rm) - vG - cross(omega,[bright;bfront;0]);
eq3= cross((omegaR3 + omegaM3 + omega),r) + cross((omegaM3 + omega),rm) - vG - cross(omega,[-bleft;-brear;0]);
eq4= cross((omegaR4 + omegaM4 + omega),r) + cross((omegaM4 + omega),rm) - vG - cross(omega,[bright;-brear;0]);            

%%
%Calcolo vettore soluzione [vGx vGy omegaz omegaR11 omegaR22 omegaR32 omegaR41 omegaM4x]
soluzione=solve([eq1(1:2);eq2(1:2);eq3(1:2);eq4(1:2)],[vGx;vGy;omegaz;omegaR11;omegaR22;omegaR32;omegaR41;omegaM4x]);

vGloc=vG;
omegaloc=omega;

%Sostituzione simbolica
vGloc = subs(vGloc,[vGx vGy],[soluzione.vGx soluzione.vGy]);
omegaloc = subs(omegaloc,omegaz,soluzione.omegaz);

%%
%Sostituzione simbolica e calcolo leggi orarie
variabili_ = [Rt,Rr,Rm,a14,a23,bfront,brear,bleft,bright];
valori_ = [0.08 0.019 0.061 pi/4 3*pi/4 0.5 0.5 0.25 0.25];
omegaIn=[omegaM1x(t) omegaM2x(t) omegaM3x(t)];

%%Definizione delle velocità angolari in input
valori_omegaIn=[2*sin(t) t t];

%Calcolo angolo di rotazione
thetaInerz(t) = int(subs(subs(omegaloc(3),omegaIn,valori_omegaIn),variabili_,valori_),t);

%Calcolo vettore velocita centro geometrico nel sdr inerziale
vGloc_sub = subs(subs(vGloc,omegaIn,valori_omegaIn),variabili_,valori_);
vGinerz = Rot(thetaInerz)*vGloc_sub;

%Calcolo leggi orarie lungo X e Y
X(t)=int(vGinerz(1),t);
Y(t)=int(vGinerz(2),t);

%%
%Plotting

% Plot X(t)
figure (1)
subplot(3,1,1)
fplot(X(t),[0 10]);
xlabel('t[s]');
ylabel('X (G) [m]');

% Plot Y(t)
subplot(3,1,2)
fplot(Y(t),[0 10]);
xlabel('t[s]');
ylabel('Y (G) [m]');

%Plot Theta(t)
subplot(3,1,3)
fplot(thetaInerz(t),[0 10]);
xlabel('t[s]');
ylabel('Theta [rad]');

%Plot Traiettoria
figure (2)
fplot(X(t),Y(t),[0 10]);
xlabel('X (G) [m]');
ylabel('Y (G) [m]');

%%
%Determinazione accelerazione del centro geometrico nel sistema di
%riferimento inerziale

aGloc = diff(vGloc,t);
axGloc(t) = aGloc(1);
ayGloc(t) = aGloc(2);
aGinerz = aGloc + cross(omegaloc,vGloc);

%%
%Determinazione accelerazione angolare del veicolo nel sistema di
%riferimento inerziale

omegaPunto = diff(omegaloc,t);
omegaPuntoInerz(t)=omegaPunto(3);

%%
%Dinamica

%%
%Definizione vettori coppia motrice
Cm1 = [Cm1x, Cm1y, Cm1z];
Cm2 = [Cm2x, Cm2y, Cm2z];
Cm3 = [Cm3x, Cm3y, Cm3z];
Cm4 = [Cm4x, Cm4y, Cm4z];

%%
%Equazione di equilibrio alla rotazione per ogni singola ruota
eq1 = Cm1 - cross(k*Rt, F1*e1);
eq2 = Cm2 - cross(k*Rt, F2*e2);
eq3 = Cm3 - cross(k*Rt, F3*e2);
eq4 = Cm4 - cross(k*Rt, F4*e1);

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
%Calcolo accelerazione,accelerazione angolare del centro geometrico
%nel sistema di riferimento locale
 
acc=(F1+F2+F3+F4)/M;
accang = ( cross([-bleft;bfront;0],F1) + cross([bright;bfront;0],F2) + cross([-bleft;-brear;0],F3) + cross([bright;-brear;0],F4) )/Iz;
%%
%Calcolo 4 set di soluzioni A B C D

A = solve(subs([acc(1)-axGloc(t);acc(2)-ayGloc(t);accang(3)-omegaPuntoInerz(t)],Cm4x,0),[Cm1x,Cm2x,Cm3x]);
B = solve(subs([acc(1)-axGloc(t);acc(2)-ayGloc(t);accang(3)-omegaPuntoInerz(t)],Cm3x,0),[Cm1x,Cm2x,Cm4x]);
C = solve(subs([acc(1)-axGloc(t);acc(2)-ayGloc(t);accang(3)-omegaPuntoInerz(t)],Cm2x,0),[Cm1x,Cm3x,Cm4x]);
D = solve(subs([acc(1)-axGloc(t);acc(2)-ayGloc(t);accang(3)-omegaPuntoInerz(t)],Cm1x,0),[Cm2x,Cm3x,Cm4x]);

%Definizione valori da sostituire nelle variabili simboliche
variabili = [Rt Rr a14 a23 M Iz bfront brear bleft bright];
valori = [0.076 0.019 pi/4 3*pi/4 141+360 26.0938 0.375 0.375 0.125 0.125];

%%
%Caso A ----- F1=0
Cm1x_(t) = 0*t;
Cm2x_(t) = subs(subs(D.Cm2x,omegaIn,valori_omegaIn),variabili,valori);
Cm3x_(t) = subs(subs(D.Cm3x,omegaIn,valori_omegaIn),variabili,valori);
Cm4x_(t) = subs(subs(D.Cm4x,omegaIn,valori_omegaIn),variabili,valori);


figure(3)
subplot(2,1,1)
title('Caso A : F1->0')
hold on
box on
fplot(Cm1x_(t),[0 10])
fplot(Cm2x_(t),[0 10])
fplot(Cm3x_(t),[0 10])
fplot(Cm4x_(t),[0 10])
legend('C_1','C_2','C_3','C_4')
xlabel('t(s)')
ylabel('coppia(Nm)')

%%
%Caso B ----- F2=0
Cm1x_(t) = subs(subs(C.Cm1x,omegaIn,valori_omegaIn),variabili,valori);
Cm2x_(t) = 0*t;
Cm3x_(t) = subs(subs(C.Cm3x,omegaIn,valori_omegaIn),variabili,valori);
Cm4x_(t) = subs(subs(C.Cm4x,omegaIn,valori_omegaIn),variabili,valori);


subplot(2,1,2)
title('Caso B : F2->0')
hold on
box on
fplot(Cm1x_(t),[0 10])
fplot(Cm2x_(t),[0 10])
fplot(Cm3x_(t),[0 10])
fplot(Cm4x_(t),[0 10])
legend('C_1','C_2','C_3','C_4')
xlabel('t(s)')
ylabel('coppia(Nm)')

%%
%Caso C ----- F3=0
Cm1x_(t) = subs(subs(B.Cm1x,omegaIn,valori_omegaIn),variabili,valori);
Cm2x_(t) = subs(subs(B.Cm2x,omegaIn,valori_omegaIn),variabili,valori);
Cm3x_(t) = 0*t;
Cm4x_(t) = subs(subs(B.Cm4x,omegaIn,valori_omegaIn),variabili,valori);

figure(4)
subplot(2,1,1)
title('Caso C : F3->0')
hold on
box on
fplot(Cm1x_(t),[0 10])
fplot(Cm2x_(t),[0 10])
fplot(Cm3x_(t),[0 10])
fplot(Cm4x_(t),[0 10])
legend('C_1','C_2','C_3','C_4')
xlabel('t(s)')
ylabel('coppia(Nm)')

%%
%Caso D ----- F4=0
Cm1x_(t) = subs(subs(A.Cm1x,omegaIn,valori_omegaIn),variabili,valori);
Cm2x_(t) = subs(subs(A.Cm2x,omegaIn,valori_omegaIn),variabili,valori);
Cm3x_(t) = subs(subs(A.Cm3x,omegaIn,valori_omegaIn),variabili,valori);
Cm4x_(t) = 0*t;

subplot(2,1,2)
title('Caso D : F4->0')
hold on
box on
fplot(Cm1x_(t),[0 10])
fplot(Cm2x_(t),[0 10])
fplot(Cm3x_(t),[0 10])
fplot(Cm4x_(t),[0 10])
legend('C_1','C_2','C_3','C_4')
xlabel('t(s)')
ylabel('coppia(Nm)')























            