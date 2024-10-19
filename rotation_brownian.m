%%%%%%%%%%%%%%%%%%%%% 1 agglomerate%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% second order stochastic non linear runge kutha %%%%%%%%%%%
%%% with dimention %%%%%

clear all;
close all;
clc

t1=cputime;

dfac=1e-6;
Kn=1e12;
skip=1000;
skip1=1000;

%constants
alpha1=0.136713;
alpha2=0.863287;
beta1=-1.512997;
beta2=1.112094;
c2=0.579182;
d2=1.18816;
a21=0.579182;
b21=-1.512997;
e21=1.18816;
g21=2.16704;
q1=0.25301;
q2=0.34026;

sph_coord_A=importdata('center_sphere.txt');
NA=size(sph_coord_A,1);

%calculation of diffusion constant
a=1e-6;
b=2e-6;
c=3e-6;
dens=1000;
mu=1.82e-5;
A_1r=3.930;
A_2r=2.058;
A_3r=0.3277;
Kb=1.38e-23;
Temp=300;
Cr_Kn=1+Kn.*(A_1r+A_2r*exp(-(A_3r/Kn)));
num=8*pi*mu*(a.^3);
den_1=((1+5.988*Kn)./Cr_Kn);
den_2=(1./((0.713*(NA.^1.63))+0.287))+((5.988*Kn).*(1./((1.184*(NA.^2.02))-0.184)));
den=den_1.*(1./(den_2));
f_r=num*den;
d_r=(Kb*Temp)./f_r;


%assigning positions and masses of particles
for i=1:1:NA   
    
    XA(i)=sph_coord_A(i,1)*a;
    YA(i)=sph_coord_A(i,2)*a;
    ZA(i)=sph_coord_A(i,3)*a;
    MA(i)=1e-15;
end

%calculating the COM of A
COM_XA=0;
COM_YA=0;
COM_ZA=0;
MA_TOTAL=0;
              
for i=1:1:NA
    COM_XA=COM_XA+XA(i)*MA(i);
    COM_YA=COM_YA+YA(i)*MA(i);
    COM_ZA=COM_ZA+ZA(i)*MA(i);
    MA_TOTAL=MA_TOTAL+MA(i); 
end
            
    COM_XA=COM_XA/MA_TOTAL;
    COM_YA=COM_YA/MA_TOTAL;
    COM_ZA=COM_ZA/MA_TOTAL;
    
for i=1:1:NA
    REL_XA(i)=XA(i)-COM_XA;
    REL_YA(i)=YA(i)-COM_YA;
    REL_ZA(i)=ZA(i)-COM_ZA;
end
for i=1:NA
    RELA{i}=[REL_XA(i);REL_YA(i);REL_ZA(i)];     %matrix of relative distances in x,y,z direction
end

%%% calculation of MOI of A  %%%
IA=[(1/5)*MA(1)*(b^2+c^2);(1/5)*MA(1)*(c^2+a^2);(1/5)*MA(1)*(a^2+b^2)];

diff_cons=((Kb*Temp)/((IA(1).^2+IA(2).^2+IA(3).^2)^0.5))^0.5;
time_period=1/diff_cons;

%initial angular velocity
omega_1=(((Kb*Temp)./((IA(1)^2)+(IA(2)^2)+(IA(3)^2))^0.5)^0.5)*randn;
omega_2=(((Kb*Temp)./((IA(1)^2)+(IA(2)^2)+(IA(3)^2))^0.5)^0.5)*randn;
omega_3=(((Kb*Temp)./((IA(1)^2)+(IA(2)^2)+(IA(3)^2))^0.5)^0.5)*randn;
omegaA=[omega_1;omega_2;omega_3;0];


%initial torque T=-f_r*omega
torque_1=-f_r*omega_1;
torque_2=-f_r*omega_2;
torque_3=-f_r*omega_3;
T=[torque_1;torque_2;torque_3];

%calculating initial angular velocity and kinetic energy
ang_vel=(omegaA(1).^2+omegaA(2).^2+omegaA(3).^2)^0.5;
ke=0.5*(IA(1)*(omegaA(1).^2)+IA(2)*(omegaA(2).^2)+IA(3)*(omegaA(3).^2));

% Initializing Euler angle in degree
r4=rand(6,1);
thetaA=180*r4(1,:);
phiA=180*r4(2,:);
psieA=180*r4(3,:);

thetaA_0=thetaA;
phiA_0=phiA;
psieA_0=psieA;

tthetaA(1,:)=thetaA;
pphiA(1,:)=phiA;
ppsieA(1,:)=psieA;

%Euler angle in radian
thetaA=thetaA*(pi/180);
phiA=phiA*(pi/180);
psieA=psieA*(pi/180);

%calculating quaternions
q1A=sin(thetaA/2)*cos((phiA-psieA)/2);
q2A=sin(thetaA/2)*sin((phiA-psieA)/2);
q3A=cos(thetaA/2)*sin((phiA+psieA)/2);
q4A=cos(thetaA/2)*cos((phiA+psieA)/2);

%differential equations
%for q1
f1A = @(t,x1A,y1A,T) y1A;
g1A = @(t,x1A,y1A,T) 0;
f2A = @(t,x1A,y1A,x2A,y2A,x3A,y3A,x4A,y4A,T) 0.5.*((x4A.*(1/IA(1))).*(T(1)+4.*(IA(2)-IA(3)).*(-y1A.*x3A+y2A.*x4A+y3A.*x1A-y4A.*x2A).*(y1A.*x2A-y2A.*x1A+y3A.*x4A-y4A.*x3A))-(x3A.*(1./IA(2))).*(T(2)+4.*(IA(3)-IA(1)).*(+y1A.*x2A-y2A.*x1A+y3A.*x4A-y4A.*x3A).*(y1A.*x4A+y2A.*x3A-y3A.*x2A-y4A.*x1A))+(x2A.*(1./IA(3))).*(T(3)+4.*(IA(1)-IA(2)).*(y1A.*x4A+y2A.*x3A-y3A.*x2A-y4A.*x1A).*(-y1A.*x3A+y2A.*x4A+y3A.*x1A-y4A.*x2A))+x1A.*((-2).*((y1A.^2)+(y2A.^2)+(y3A.^2)+(y4A.^2))));
%for q2
f3A = @(t,x2A,y2A,T) y2A;
g3A = @(t,x2A,y2A,T) 0;
f4A = @(t,x1A,y1A,x2A,y2A,x3A,y3A,x4A,y4A,T) 0.5.*(x3A.*(1/IA(1)).*(T(1)+4*(IA(2)-IA(3)).*(-y1A.*x3A+y2A.*x4A+y3A.*x1A-y4A.*x2A).*(y1A.*x2A-y2A.*x1A+y3A.*x4A-y4A.*x3A))+x4A.*(1./IA(2)).*(T(2)+4*(IA(3)-IA(1)).*(+y1A.*x2A-y2A.*x1A+y3A.*x4A-y4A.*x3A).*(y1A.*x4A+y2A.*x3A-y3A.*x2A-y4A.*x1A))-x1A.*(1./IA(3)).*(T(3)+4.*(IA(1)-IA(2)).*(y1A.*x4A+y2A.*x3A-y3A.*x2A-y4A.*x1A).*(-y1A.*x3A+y2A.*x4A+y3A.*x1A-y4A.*x2A))+x2A.*((-2).*((y1A.^2)+(y2A.^2)+(y3A.^2)+(y4A.^2))));
%for q3
f5A = @(t,x3A,y3A,T) y3A;
g5A = @(t,x3A,y3A,T) 0;
f6A = @(t,x1A,y1A,x2A,y2A,x3A,y3A,x4A,y4A,T) 0.5.*(-x2A.*(1/IA(1)).*(T(1)+4*(IA(2)-IA(3)).*(-y1A.*x3A+y2A.*x4A+y3A.*x1A-y4A.*x2A).*(y1A.*x2A-y2A.*x1A+y3A.*x4A-y4A.*x3A))+x1A.*(1./IA(2)).*(T(2)+4*(IA(3)-IA(1)).*(+y1A.*x2A-y2A.*x1A+y3A.*x4A-y4A.*x3A).*(y1A.*x4A+y2A.*x3A-y3A.*x2A-y4A.*x1A))+x4A.*(1./IA(3)).*(T(3)+4.*(IA(1)-IA(2)).*(y1A.*x4A+y2A.*x3A-y3A.*x2A-y4A.*x1A).*(-y1A.*x3A+y2A.*x4A+y3A.*x1A-y4A.*x2A))+x3A.*((-2).*((y1A.^2)+(y2A.^2)+(y3A.^2)+(y4A.^2))));
%for q4
f7A = @(t,x4A,y4A,T) y4A;
g7A = @(t,x4A,y4A,T) 0;
f8A = @(t,x1A,y1A,x2A,y2A,x3A,y3A,x4A,y4A,T) 0.5.*(-x1A.*(1/IA(1)).*(T(1)+4*(IA(2)-IA(3)).*(-y1A.*x3A+y2A.*x4A+y3A.*x1A-y4A.*x2A).*(y1A.*x2A-y2A.*x1A+y3A.*x4A-y4A.*x3A))-x2A.*(1./IA(2)).*(T(2)+4*(IA(3)-IA(1)).*(+y1A.*x2A-y2A.*x1A+y3A.*x4A-y4A.*x3A).*(y1A.*x4A+y2A.*x3A-y3A.*x2A-y4A.*x1A))-x3A.*(1./IA(3)).*(T(3)+4.*(IA(1)-IA(2)).*(y1A.*x4A+y2A.*x3A-y3A.*x2A-y4A.*x1A).*(-y1A.*x3A+y2A.*x4A+y3A.*x1A-y4A.*x2A))+x4A.*((-2).*((y1A.^2)+(y2A.^2)+(y3A.^2)+(y4A.^2))));

%relation between omega nd q, q_dot
WA=[q4A q3A -q2A -q1A;-q3A q4A q1A -q2A;q2A -q1A q4A -q3A;q1A q2A q3A q4A] ;
WA_inv=inv(WA);
qA_dot=0.5.*(WA_inv*omegaA);
dA=1;
dA_t=1.00000;

%initial value of q and qdot
x1A=q1A;
y1A=qA_dot(1);
x2A=q2A;
y2A=qA_dot(2);
x3A=q3A;
y3A=qA_dot(3);
x4A=q4A;
y4A=qA_dot(4);

xx1A(1,:)=x1A;
yy1A(1,:)=y1A;
xx2A(1,:)=x2A;
yy2A(1,:)=y2A;
xx3A(1,:)=x3A;
yy3A(1,:)=y3A;
xx4A(1,:)=x4A;
yy4A(1,:)=y4A;

v= x1A*y1A+x2A*y2A+x3A*y3A+x4A*y4A;
qsq=(x1A.^2)+(x2A.^2)+(x3A.^2)+(x4A.^2);
 
vv(1,:)= v;
qqsq(1,:)=qsq;
ddA(1,:)=dA;
oomega_1(1,:)=omega_1;
oomega_2(1,:)=omega_2;
oomega_3(1,:)=omega_3;
aang_vel(1,:)=ang_vel;
kke(1,:)=ke;
    
oomega_1(1)=omega_1;
oomega_2(1)=omega_2;
oomega_3(1)=omega_3;
aang_vel(1)=ang_vel;
kke(1)=ke;
xx_1A(1)=x1A;
xx_2A(1)=x2A;
xx_3A(1)=x3A;
xx_4A(1)=x4A;
xx1A_cor(1)=x1A;
xx2A_cor(1)=x2A;
xx3A_cor(1)=x3A;
xx4A_cor(1)=x4A;

n(1)=0;
nn=1;
nnn=1;
i=1;
j=1;
t=0;
tt(1)=0;
time=t(1);
ttr=((IA(1).^2+(IA(2).^2)+(IA(3)).^2)^0.5)./f_r;

while t<=5*time_period

%calculating timestep
h1=1/y1A;
h2=1/y2A;
h3=1/y3A;
h4=1/y4A;
h5=1/(abs(f2A(t,x1A,y1A,x2A,y2A,x3A,y3A,x4A,y4A,T)))^0.5;
h6=1/(abs(f4A(t,x1A,y1A,x2A,y2A,x3A,y3A,x4A,y4A,T)))^0.5;
h7=1/(abs(f6A(t,x1A,y1A,x2A,y2A,x3A,y3A,x4A,y4A,T)))^0.5;
h8=1/(abs(f8A(t,x1A,y1A,x2A,y2A,x3A,y3A,x4A,y4A,T)))^0.5;
h9=((IA(1)^2)/(f_r*Kb*Temp))^(1/3);
h10=((IA(2)^2)/(f_r*Kb*Temp))^(1/3);
h11=((IA(3)^2)/(f_r*Kb*Temp))^(1/3);
hhh=[h1;h2;h3;h4;h5;h6;h7;h8;h9;h10;h11];
hhh_abs(:,1) = abs(hhh(:,1));
hhhh=min(hhh_abs);
h=dfac*hhhh;

%generating random number for stochastic components of torque
sd1=(q1*f_r*Kb*Temp*h/(IA(1)^2))^0.5;
sd2=(q1*f_r*Kb*Temp*h/(IA(2)^2))^0.5;
sd3=(q1*f_r*Kb*Temp*h/(IA(3)^2))^0.5;
sd11=(q2*f_r*Kb*Temp*h/(IA(1)^2))^0.5;
sd22=(q2*f_r*Kb*Temp*h/(IA(2)^2))^0.5;
sd33=(q2*f_r*Kb*Temp*h/(IA(3)^2))^0.5;

w1=normrnd(0,sd1)*100;
w2=normrnd(0,sd2)*100;
w3=normrnd(0,sd3)*100;
w11=normrnd(0,sd11)*100;
w22=normrnd(0,sd22)*100;
w33=normrnd(0,sd33)*100;

  %first order
  % for A
  %for q1
  k1x1A = h*f1A(t,x1A,y1A,T);
  j1x1A=0;
  k1y1A = h*f2A(t,x1A,y1A,x2A,y2A,x3A,y3A,x4A,y4A,T);
  j1y1A = 0.5*(x4A*w1-x3A*w2+x2A*w3);
  %for q2
  k1x2A = h*f3A(t,x2A,y2A,T);
  j1x2A=0;
  k1y2A = h*f4A(t,x1A,y1A,x2A,y2A,x3A,y3A,x4A,y4A,T);
  j1y2A = 0.5*(x3A*w1+x4A*w2-x1A*w3);
  %for q3
  k1x3A = h*f5A(t,x3A,y3A,T);
  j1x3A=0;
  k1y3A = h*f6A(t,x1A,y1A,x2A,y2A,x3A,y3A,x4A,y4A,T);
  j1y3A = 0.5*(-x2A*w1+x1A*w2+x4A*w3);
  %for q4
  k1x4A = h*f7A(t,x4A,y4A,T);
  j1x4A=0;
  k1y4A = h*f8A(t,x1A,y1A,x2A,y2A,x3A,y3A,x4A,y4A,T);
  j1y4A = 0.5*(-x1A*w1-x2A*w2-x3A*w3);

  %second order
  % for A
  %for q1
  k2x1A = h*f1A(t+c2*h, x1A+a21*k1x1A+b21*j1x1A, y1A,T);
  j2x1A=0;
  k2y1A = h*f2A(t+c2*h, x1A, y1A+a21*k1y1A+b21*j1y1A, x2A, y2A, x3A, y3A, x4A, y4A, T);
  j2y1A = 0.5*(x4A*w11-x3A*w22+x2A*w33);
  %for q2
  k2x2A = h*f3A(t+c2*h, x2A+a21*k1x2A+b21*j1x2A, y2A,T);
  j2x2A=0;
  k2y2A = h*f4A(t+c2*h, x1A, y1A, x2A, y2A+a21*k1y2A+b21*j1y2A, x3A, y3A, x4A, y4A, T);
  j2y2A = 0.5*(x3A*w11+x4A*w22-x1A*w33);
  %for q3
  k2x3A = h*f5A(t+c2*h, x3A+a21*k1x3A+b21*j1x3A, y3A,T);
  j2x3A=0;
  k2y3A = h*f6A(t+c2*h, x1A, y1A, x2A, y2A, x3A, y3A+a21*k1y3A+b21*j1y3A, x4A, y4A, T);  
  j2y3A = 0.5*(-x2A*w11+x1A*w22+x4A*w33);
  %for q4
  k2x4A = h*f7A(t+c2*h, x4A+a21*k1x4A+b21*j1x4A, y4A,T);
  j2x4A=0;
  k2y4A = h*f8A(t+c2*h, x1A, y1A, x2A, y2A, x3A, y3A, x4A, y4A+a21*k1y4A+b21*j1y4A, T);
  j2y4A = 0.5*(-x1A*w11-x2A*w22-x3A*w33);
  

   %next iteration values
  %q1
  x1A = x1A+alpha1*k1x1A+beta1*j1x1A+alpha2*k2x1A+beta2*j2x1A;
  y1A = y1A+alpha1*k1y1A+beta1*j1y1A+alpha2*k2y1A+beta2*j2y1A;

  %q2
  x2A = x2A+alpha1*k1x2A+beta1*j1x2A+alpha2*k2x2A+beta2*j2x2A;
  y2A = y2A+alpha1*k1y2A+beta1*j1y2A+alpha2*k2y2A+beta2*j2y2A;

  %q3
  x3A = x3A+alpha1*k1x3A+beta1*j1x3A+alpha2*k2x3A+beta2*j2x3A;
  y3A = y3A+alpha1*k1y3A+beta1*j1y3A+alpha2*k2y3A+beta2*j2y3A;

  %q4
  x4A = x4A+alpha1*k1x4A+beta1*j1x4A+alpha2*k2x4A+beta2*j2x4A;
  y4A = y4A+alpha1*k1y4A+beta1*j1y4A+alpha2*k2y4A+beta2*j2y4A;


v= x1A*y1A+x2A*y2A+x3A*y3A+x4A*y4A;
ymat=[y1A;y2A;y3A;y4A];
WA=[x4A x3A -x2A -x1A;-x3A x4A x1A -x2A;x2A -x1A x4A -x3A;x1A x2A x3A x4A] ;
omegaA=2.*(WA*ymat);

%updating angular velocity
ang_vel=((omegaA(1)).^2+(omegaA(2)).^2+(omegaA(3)).^2)^0.5;
omega_1= omegaA(1);
omega_2= omegaA(2);
omega_3= omegaA(3);

%updating kinetic energy
ke=0.5*(IA(1)*(omegaA(1).^2)+IA(2)*(omegaA(2).^2)+IA(3)*(omegaA(3).^2));

%updating kinetic energy
torque_1=-f_r*omegaA(1);
torque_2=-f_r*omegaA(2);
torque_3=-f_r*omegaA(3);

T=[torque_1;torque_2;torque_3];

%updating rotation matrix
 R11A=(x1A.^2)+(x4A.^2)-0.5;
R12A=x1A.*x2A+x3A.*x4A;
R13A=x1A.*x3A-x2A.*x4A;
R21A=x1A.*x2A-x3A.*x4A;
R22A=(x2A.^2)+(x4A.^2)-0.5;
R23A=x2A.*x3A+x1A.*x4A;
R31A=x1A.*x3A+x2A.*x4A;
R32A=x2A.*x3A-x1A.*x4A;
R33A=(x3A.^2)+(x4A.^2)-0.5;
  
RA=2.*[R11A R12A R13A;R21A R22A R23A;R31A R32A R33A];
dA=det(RA);
dA=double(dA);  

sin_theta=2*((((x1A^2)+(x2A^2))*(1-(x1A^2)-(x2A^2)))^0.5);
cos_theta=1-2*((x1A^2)+(x2A^2));
sin_phi=2*(x1A*x3A+x2A*x4A)/sin_theta;
cos_phi=2*(x1A*x4A-x2A*x3A)/sin_theta;
sin_psi=2*(x1A*x3A-x2A*x4A)/sin_theta;
cos_psi=2*(x1A*x4A+x2A*x3A)/sin_theta;

thetaA=acos(cos_theta); %in radian
phiA=acos(cos_phi);
psieA=acos(cos_psi);

thetaA=thetaA*(180/pi);  %in degree
phiA=phiA*(180/pi);
psieA=psieA*(180/pi);


  if mod(i+1,skip)==0 %printing values
      nn=nn+1;
     ddA(nn,:)=dA;
     RRA{nn}=RA;
     h_s(nn,:)=h;
     writematrix(ddA,'determinant.txt');
     writematrix(h_s,'timestep.txt');     

    oomega_1(nn,:)=omega_1;
    oomega_2(nn,:)=omega_2;
    oomega_3(nn,:)=omega_3;
    aang_vel(nn,:)=ang_vel;

    writematrix(oomega_1,'ang_vel_x.txt');
    writematrix(oomega_2,'ang_vel_y.txt');
    writematrix(oomega_3,'ang_vel_z.txt');
    writematrix(aang_vel,'mag_ang_vel.txt');

    kke(nn,:)=ke;
    writematrix(kke,'kin_energy.txt');

    xx1A(nn,:)=x1A;
    xx2A(nn,:)=x2A;
    xx3A(nn,:)=x3A;
    xx4A(nn,:)=x4A;
    yy1A(nn,:)=y1A;
    yy2A(nn,:)=y2A;
    yy3A(nn,:)=y3A;
    yy4A(nn,:)=y4A;

    tt(nn,:)=t;
    writematrix(xx1A,'q1.txt');
    writematrix(xx2A,'q2.txt');
    writematrix(xx3A,'q3.txt');
    writematrix(xx4A,'q4.txt');
    writematrix(yy1A,'q1dot.txt');
    writematrix(yy2A,'q2dot.txt');
    writematrix(yy3A,'q3dot.txt');
    writematrix(yy4A,'q4dot.txt');
    writematrix(tt,'time.txt');
   
    tthetaA(nn,:)=thetaA;
    pphiA(nn,:)=phiA;
    ppsieA(nn,:)=psieA;

    writematrix(tthetaA,'theta.txt');
    writematrix(pphiA,'phi.txt');
    writematrix(ppsieA,'psie.txt');
    
  end

  if mod(i+1,skip1)==0
       RRA{nnn}=RA;
%       cla;
     for i_sphA=1:size(sph_coord_A,1)
     tempA{i_sphA}=sph_coord_A(i_sphA,:);
     tempA{i_sphA}=RRA{nn}*tempA{i_sphA}';
    tempA_print{i_sphA}=tempA{i_sphA}';%+pos_A{i+1};
    x(nnn,:)=tempA{i_sphA}(1);
    y(nnn,:)=tempA{i_sphA}(2);
    z(nnn,:)=tempA{i_sphA}(3);

    writematrix(x,'x.txt');
    writematrix(y,'y.txt');
    writematrix(z,'z.txt');

    nnn=nnn+1;

     end
  end

    t=t+h;
    i=i+1;
    t2=cputime-t1;
    writematrix(t2,'t_cpu.txt');
    writematrix(t,'t_sim.txt');

 end
