
clc
clear all 
close all 


% Parameters of Blasius Equation ------------------------------------------

rho=1;
mu=1.8e-5;
U=0.18;
L=1;
h_x=0.01;
eta_i=0;
eta_f=7 ;
nu=mu/rho;

% initial deta for Blasius Equation ---------------------------------------
h_eta=0.01;
y(1,1)=0     ;  
y(1,2)=0     ;
f_pp_f=1;
error=0.00000001 ;

eta=[eta_i:h_eta:eta_f];
n=length(eta);

% Shooting Method Algorithm for Blasius Equation --------------------------

z(1)=0;
y(1,3)=z(1) ;
y=rk4(eta,y);
ph(1)=y(end,2);

z(2)=1;
y(1,3)= z(2) ;
y=rk4(eta,y);
ph(2)=y(end,2);
for j=2:1:100
z(j+1)=z(j-1)+( f_pp_f-ph(j-1) )*( z(j)-z(j-1) )/( ph(j)-ph(j-1) ) ;  
y(1,3)=z(j+1);
y=rk4(eta,y);
ph(j+1)=y(end,2); 
if( abs(y(end,2)-f_pp_f)<error )
break   
end  
end  

% figure(1) ---------------------------------------------------------------  

figure(1)
plot(eta,y(:,1), eta, y(:,2), eta, y(:,3), 'LineWidth', 2)

xlim([0 5])
ylim([0 2])
title('Solution of Blasius eqution', 'FontSize', 12);
ylabel('f, f'' , f''''', 'FontSize', 15);
xlabel('\eta', 'FontSize', 20);
grid on
Legend1 = {'f(\eta)', 'f''(\eta)', 'f''''(\eta)'};
legend(Legend1, 'FontSize', 14);



r(:,1)=eta;
r(:,2:1:4)=y(:,1:1:3);
disp(['    eta         f         f_p      f_pp'         ])
disp(['------------------------------------------------'])      
disp(r)
disp(['________________________________________________'])
              
        
        




% delta  ------------------------------------------------------------------

x=0:h_x:L;
for i = 1:length(x)
delta(i) = 5*sqrt(x(i)*nu/U);
delta_star(i)=delta(i)/3;
theta(i)=delta_star(i)/3;
end
figure(2)
plot(x,delta, x, delta_star, x, theta, 'LineWidth', 2)  
   
ylabel('\delta, \delta^* , \theta', 'FontSize', 15);
xlabel('x', 'FontSize', 20);
Legend1 = {'\delta', '\delta^*' , '\theta'};
legend(Legend1, 'FontSize', 14);


clear r
r(:,1)=x;
r(:,2)=delta;
r(:,3)=delta_star;
r(:,4)=theta;

disp(['    x        delta    delta_star  theta'   ])
disp(['------------------------------------------------'])      
disp(r)
disp(['________________________________________________'])






% tau  -------------------------------------------------------------------

for i = 1 : length(x)
tau(i)=y(1,3)*mu*U*sqrt(U/(2*nu*x(i)));
end

figure(3)
plot(x,tau, 'LineWidth', 2)  
xlabel('x (m)')
ylabel('Shear stress (Pa)')
title(['U_i_n_f = ' num2str(U)])




clear r
r(:,1)=x;
r(:,2)=tau;
disp(['    x         tau         '   ])
disp(['------------------------------------------------'])      
disp(r)
disp(['________________________________________________'])




% u   ------------------------------------------------------------------


figure(4)

eta=eta';
loc_x = [0.4 0.8 1.2 1.6 ];

u=U*y(:,2);

clear r
r(:,1)=u;

for i=1:1:length(loc_x) 
    y1=eta.*sqrt(nu*loc_x(i)/U);    
    plot(u,y1,  'LineWidth', 2)
    hold all
    

r(:,i+1)=y1;
    
end

xlabel('u (m/s)')
ylabel('y (m)')
Legend1 = {' x = 0.4', ' x = 0.8', ' x = 1.2',' x = 1.6'};
legend(Legend1, 'FontSize', 13);




disp(['     u     y(x=0.4)   y(x=0.8)  y(x=1.2)  y(x=1.6)' ])
disp(['---------------------------------------------------'])      
disp(r)
disp(['___________________________________________________'])




% v   ------------------------------------------------------------------


figure(5)
clear r
 for i=1:1:length(loc_x) 
     y1=eta*sqrt(nu*loc_x(i)/U);
     v=0.5*sqrt(nu*U/loc_x(i))*(eta.*y(:,2)-y(:,1));
     plot(v,y1,  'LineWidth', 2)
     hold all
     
r(:,2*i-1)=v;
r(:,2*i)=y1;
   
 end

xlabel('v (m/s)')
ylabel('y (m)')
Legend1 = {' x = 0.4', ' x = 0.8', ' x = 1.2',' x = 1.6'};
legend(Legend1, 'FontSize', 13);
 
disp(['         x=0.4               x=0.8               x=1.2              x=1.6      ']) 
disp(['    _________________   ________________    ________________   ________________'])
disp(['       v          y       v         y         v         y         v         y  '])
disp(['-------------------------------------------------------------------------------'])      
disp(r)
disp(['________________________________________________________________________________'])






% function ----------------------------------------------------------------
 
function y=rk4(x,y)
n=length(x);
h=(x(end)-x(1))/(n-1);
Y(:)=y(1,:);

for i=1:1:n-1
    K1=h*func(x(i),Y);
    K2=h*func(x(i)+h/2,Y+K1/2);
    K3=h*func(x(i)+h/2,Y+K2/2);
    K4=h*func(x(i)+h,Y+K3);
    Y=Y+(K1+2*K2+2*K3+K4)/6;
    y(i+1,:)=Y;
end
function f=func(X,Y)
f(1)=Y(2);
f(2)=Y(3);
f(3)=-0.5*Y(1)*Y(3);
end

end







