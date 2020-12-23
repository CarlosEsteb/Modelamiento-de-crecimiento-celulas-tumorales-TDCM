clear; clc;
close all;
  %Tiempos
  t0 = 0;
  tf =50;
  h =0.01;
  t =t0:h:tf;
  
  %Condiciones iniciales
  x0 = [0.5,0.2,0.1];
  
  function [xdot] = nonlinear(t,x)
  global fileID

%Paramters
 k3=1;A12=1;A13=2.5;R1=0.5;
 R2=4.6;A21=1.5;R3=8;D3=0.5;A31=0.2;
 w=pi*10;
xd1 =0.1; % Xd1 deseado
xd3 =0.4; % Xd3 deseado
a1=1; 
a2=5; 
k=0.05; % --------------OJO
%erros
e2=x(3)-xd3; 
e1=x(1)-xd1;
%superficie
s=a1*(e1)+a2*(e2);
trinche=1;
ueq=(-a1/a2)*(x(1)*(1-x(1))-A12*x(1)*x(2))-((R3*x(1)*x(3))/(x(1)+k3)-A31*x(1)*x(3)-D3*x(3));
udis=-k*abs(s)*sign(s/trinche)
u=ueq+udis; %Sliding control 
%Equiations with control 
 f1=x(1)*(1-x(1))-A12*x(1)*x(2)-A13*x(1)*x(3);
 f2=R2*x(2)*(1-x(2))-A21*x(1)*x(2);
 f3=(R3*x(1)*x(3))/(x(1)+k3)-A31*x(1)*x(3)-D3*x(3)+u;
 
 
 
 xdot = [ f1;f2;f3]; %sistem control
  
  
  fprintf(fileID,'%10.6f %10.6f %10.6f %10.6f %10.6f \n',t,u,x(1),x(2),x(3))
endfunction

 function [t,xdo] = RK4(f,t0,tf,x0,h) 
  % Runge Kutta Method 4th Order 
  t = t0:h:tf; 
  length(t); 
  x=x0'; 
  xdo=x0;  
  for i=1:(length(t)-1)       
    k1 = f(t(i),x);     
    k2 = f(t(i)+0.5*h,x+0.5*h*k1);     
    k3 = f((t(i)+0.5*h),(x+0.5*h*k2));     
    k4 = f((t(i)+h),(x+k3*h));      
    x = x + (1/6)*(k1+2*k2+2*k3+k4)*h;              
    xdo(i+1,:)=x; 
  end
endfunction
  global fileID
  fileID = fopen('datos.txt','w');
  [t,xdo] = RK4(@(t,x) nonlinear(t,x),t0,tf,x0,h);
  fclose(fileID);
  hold on;
  plot(t,xdo,'linewidth',3)
  legend('Células Tumorales','Células Inmunologicas','Células Inmunologicas en reposo');
  grid minor;