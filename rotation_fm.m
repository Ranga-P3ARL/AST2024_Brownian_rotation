%%%%%%%%%%%%%%%%%%%%% 1 agglomerate%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% free mol rotation %%%%%%%%%%%

clear all
close all;
clc

t1=cputime;

dfac=1e-5;
skip=100;
skip1=200;

sph_coord_A=importdata('center_sphere.txt');
NA=size(sph_coord_A,1);

%calculation of diffusion constant
a=1e-6;
b=2e-6;
c=3e-6;
dens=1000;
Kb=1.38e-23;
Temp=300;

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
% principle moment of inertia tensor of A
IA=[(1/5)*MA(1)*(b^2+c^2);(1/5)*MA(1)*(c^2+a^2);(1/5)*MA(1)*(a^2+b^2)];

diff_cons=((Kb*Temp)/((IA(1).^2+IA(2).^2+IA(3).^2)^0.5))^0.5; 
time_period=1/diff_cons;

%initial angular velocity
omega_1=(((Kb*Temp)./((IA(1)^2)+(IA(2)^2)+(IA(3)^2))^0.5)^0.5)*randn;
omega_2=(((Kb*Temp)./((IA(1)^2)+(IA(2)^2)+(IA(3)^2))^0.5)^0.5)*randn;
omega_3=(((Kb*Temp)./((IA(1)^2)+(IA(2)^2)+(IA(3)^2))^0.5)^0.5)*randn;
omegaA=[omega_1;omega_2;omega_3;0];

%initial torque
torque_1=0;
torque_2=0;
torque_3=0;
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
f2A = @(t,x1A,y1A,x2A,y2A,x3A,y3A,x4A,y4A,T) 0.5.*(x4A.*(1/IA(1)).*(T(1)+4*(IA(2)-IA(3)).*(-y1A.*x3A+y2A.*x4A+y3A.*x1A-y4A.*x2A).*(y1A.*x2A-y2A.*x1A+y3A.*x4A-y4A.*x3A))-x3A.*(1./IA(2)).*(T(2)+4*(IA(3)-IA(1)).*(+y1A.*x2A-y2A.*x1A+y3A.*x4A-y4A.*x3A).*(y1A.*x4A+y2A.*x3A-y3A.*x2A-y4A.*x1A))+x2A.*(1./IA(3)).*(T(3)+4.*(IA(1)-IA(2)).*(y1A.*x4A+y2A.*x3A-y3A.*x2A-y4A.*x1A).*(-y1A.*x3A+y2A.*x4A+y3A.*x1A-y4A.*x2A))+x1A.*((-2).*((y1A.^2)+(y2A.^2)+(y3A.^2)+(y4A.^2))));
%for q2
f3A = @(t,x2A,y2A,T) y2A;
f4A = @(t,x1A,y1A,x2A,y2A,x3A,y3A,x4A,y4A,T) 0.5.*(x3A.*(1/IA(1)).*(T(1)+4*(IA(2)-IA(3)).*(-y1A.*x3A+y2A.*x4A+y3A.*x1A-y4A.*x2A).*(y1A.*x2A-y2A.*x1A+y3A.*x4A-y4A.*x3A))+x4A.*(1./IA(2)).*(T(2)+4*(IA(3)-IA(1)).*(+y1A.*x2A-y2A.*x1A+y3A.*x4A-y4A.*x3A).*(y1A.*x4A+y2A.*x3A-y3A.*x2A-y4A.*x1A))-x1A.*(1./IA(3)).*(T(3)+4.*(IA(1)-IA(2)).*(y1A.*x4A+y2A.*x3A-y3A.*x2A-y4A.*x1A).*(-y1A.*x3A+y2A.*x4A+y3A.*x1A-y4A.*x2A))+x2A.*((-2).*((y1A.^2)+(y2A.^2)+(y3A.^2)+(y4A.^2))));
%for q3
f5A = @(t,x3A,y3A,T) y3A;
f6A = @(t,x1A,y1A,x2A,y2A,x3A,y3A,x4A,y4A,T) 0.5.*(-x2A.*(1/IA(1)).*(T(1)+4*(IA(2)-IA(3)).*(-y1A.*x3A+y2A.*x4A+y3A.*x1A-y4A.*x2A).*(y1A.*x2A-y2A.*x1A+y3A.*x4A-y4A.*x3A))+x1A.*(1./IA(2)).*(T(2)+4*(IA(3)-IA(1)).*(+y1A.*x2A-y2A.*x1A+y3A.*x4A-y4A.*x3A).*(y1A.*x4A+y2A.*x3A-y3A.*x2A-y4A.*x1A))+x4A.*(1./IA(3)).*(T(3)+4.*(IA(1)-IA(2)).*(y1A.*x4A+y2A.*x3A-y3A.*x2A-y4A.*x1A).*(-y1A.*x3A+y2A.*x4A+y3A.*x1A-y4A.*x2A))+x3A.*((-2).*((y1A.^2)+(y2A.^2)+(y3A.^2)+(y4A.^2))));
%for q4
f7A = @(t,x4A,y4A,T) y4A;
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


n(1)=0;
nn=1;
nnn=1;
i=1;
j=1;
t=0;
tt(1,:)=0;
time=t(1);

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
hhh=[h1;h2;h3;h4;h5;h6;h7;h8];
hhh_abs(:,1) = abs(hhh(:,1));
hhhh=min(hhh_abs);
h=dfac*hhhh;

  %first order
  % for A
  k1x1A = f1A(t,x1A,y1A,T);
  k1y1A = f2A(t,x1A,y1A,x2A,y2A,x3A,y3A,x4A,y4A,T);
  
  k1x2A = f3A(t,x2A,y2A,T);
  k1y2A = f4A(t,x1A,y1A,x2A,y2A,x3A,y3A,x4A,y4A,T);
  
  k1x3A = f5A(t,x3A,y3A,T);
  k1y3A = f6A(t,x1A,y1A,x2A,y2A,x3A,y3A,x4A,y4A,T);
  
  k1x4A = f7A(t,x4A,y4A,T);
  k1y4A = f8A(t,x1A,y1A,x2A,y2A,x3A,y3A,x4A,y4A,T);
  

  %second order
  % for A
  k2x1A = f1A(t+0.5*h, x1A+0.5*k1x1A*h, y1A+0.5*k1y1A*h,T);
  k2y1A = f2A(t+0.5*h, x1A+0.5*k1x1A*h, y1A+0.5*k1y1A*h, x2A, y2A, x3A, y3A, x4A, y4A, T);
  
  k2x2A = f3A(t+0.5*h, x2A+0.5*k1x2A*h, y2A+0.5*k1y2A*h,T);
  k2y2A = f4A(t+0.5*h, x1A, y1A, x2A+0.5*k1x2A*h, y2A+0.5*k1y2A*h, x3A, y3A, x4A, y4A, T);
  
  k2x3A = f5A(t+0.5*h, x3A+0.5*k1x3A*h, y3A+0.5*k1y3A*h,T);
  k2y3A = f6A(t+0.5*h, x1A, y1A, x2A, y2A, x3A+0.5*k1x3A*h, y3A+0.5*k1y3A*h, x4A, y4A, T);
 
  k2x4A = f7A(t+0.5*h, x4A+0.5*k1x4A*h, y4A+0.5*k1y4A*h,T);
  k2y4A = f8A(t+0.5*h, x1A, y1A, x2A, y2A, x3A, y3A, x4A+0.5*k1x4A*h, y4A+0.5*k1y4A*h, T);

  %third order
  % for A
  k3x1A = f1A(t+0.5*h, x1A+0.5*k2x1A*h, y1A+0.5*k2y1A*h,T);
  k3y1A = f2A(t+0.5*h, x1A+0.5*k2x1A*h, y1A+0.5*k2y1A*h, x2A, y2A, x3A, y3A, x4A, y4A, T);
  
  k3x2A = f3A(t+0.5*h, x2A+0.5*k2x2A*h, y2A+0.5*k2y2A*h,T);
  k3y2A = f4A(t+0.5*h, x1A, y1A, x2A+0.5*k2x2A*h, y2A+0.5*k2y2A*h, x3A, y3A, x4A, y4A, T);
  
  k3x3A = f5A(t+0.5*h, x3A+0.5*k2x3A*h, y3A+0.5*k2y3A*h,T);
  k3y3A = f6A(t+0.5*h, x1A, y1A, x2A, y2A, x3A+0.5*k2x3A*h, y3A+0.5*k2y3A*h, x4A, y4A, T);
  
  k3x4A = f7A(t+0.5*h, x4A+0.5*k2x4A*h, y4A+0.5*k2y4A*h,T);
  k3y4A = f8A(t+0.5*h, x1A, y1A, x2A, y2A, x3A, y3A, x4A+0.5*k2x4A*h, y4A+0.5*k2y4A*h, T);

  %forth order
  % for A
  k4x1A = f1A(t+h, x1A+k3x1A*h, y1A+k3y1A*h, T);
  k4y1A = f2A(t+h, x1A+k3x1A*h, y1A+k3y1A*h, x2A, y2A, x3A, y3A, x4A, y4A, T);
  
  k4x2A = f3A(t+h, x2A+k3x2A*h, y2A+k3y2A*h, T);
  k4y2A = f4A(t+h, x1A, y1A, x2A+k3x2A*h, y2A+k3y2A*h, x3A, y3A, x4A, y4A, T);
  
  k4x3A = f5A(t+h, x3A+k3x3A*h, y3A+k3y3A*h, T);
  k4y3A = f6A(t+h, x1A, y1A, x2A, y2A, x3A+k3x3A*h, y3A+k3y3A*h, x4A, y4A, T);
  
  k4x4A = f7A(t+h, x4A+k3x4A*h, y4A+k3y4A*h,T);
  k4y4A = f8A(t+h, x1A, y1A, x2A, y2A, x3A, y3A, x4A+k3x4A*h, y4A+k3y4A*h, T);
  

  %next iteration values
%for A
  x1A = x1A+1/6*(k1x1A+2*k2x1A+2*k3x1A+k4x1A)*h;
  y1A = y1A+1/6*(k1y1A+2*k2y1A+2*k3y1A+k4y1A)*h;
  
  x2A = x2A+1/6*(k1x2A+2*k2x2A+2*k3x2A+k4x2A)*h;
  y2A = y2A+1/6*(k1y2A+2*k2y2A+2*k3y2A+k4y2A)*h;
  
  x3A = x3A+1/6*(k1x3A+2*k2x3A+2*k3x3A+k4x3A)*h;
  y3A = y3A+1/6*(k1y3A+2*k2y3A+2*k3y3A+k4y3A)*h;
  
  x4A = x4A+1/6*(k1x4A+2*k2x4A+2*k3x4A+k4x4A)*h;
  y4A = y4A+1/6*(k1y4A+2*k2y4A+2*k3y4A+k4y4A)*h;


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

%updating torque
torque_1=0;
torque_2=0;
torque_3=0;

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

%updating Euler angles
sin_theta=2*((((x1A^2)+(x2A^2))*(1-(x1A^2)-(x2A^2)))^0.5);
cos_theta=1-2*((x1A^2)+(x2A^2));
sin_phi=2*(x1A*x3A+x2A*x4A)/sin_theta;
cos_phi=2*(x1A*x4A-x2A*x3A)/sin_theta;
sin_psi=2*(x1A*x3A-x2A*x4A)/sin_theta;
cos_psi=2*(x1A*x4A+x2A*x3A)/sin_theta;

thetaA=acos(cos_theta); %in radian
phiA=acos(cos_phi);
psieA=acos(cos_psi);

thetaA=thetaA*(180/pi); %in degree
phiA=phiA*(180/pi);
psieA=psieA*(180/pi);

  if mod(i+1,skip)==0  %printing values
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
    tempA_print{i_sphA}=tempA{i_sphA}';
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

