%   ****************************************************************************************************************%
%                       %% BRACKET OPERATOR PENALTY METHOD  - PHASE III %%
%                       %% POWELL CONJUGATE DIRECTION METHOD - PHASE II %%
%                       %% BOUNDING PHASE AND SECANT METHOD  - PHASE I  %%
%                       %% GROUP (G-15) : AVINASH CHETRY 206103107 , SASANKA SEKHAR SAIKIA 214103429
%   ****************************************************************************************************************%
close all;
clc;
clear all;
format short;
%% USER DEFINED VARIABLES
N = input("Enter the dimension of Design Vector:\n");
%% USER DEFINED INPUT VECTOR
x0=zeros(N,1);
for j=1:N
    x0(j)=input('Enter the Elements: ');  % coordinates of the first attempt (column vector)
end
%disp(x0);
%% MAIN CODE Penalty Function Method Phase 3
tol1 = 10e-5;       % tolerance or accuracy
tol2 = 10e-5;
R0 = 0.1;           % initial Penalty parameter
% x0 = [0 0];
c = input('enter 0 for interior penalty, 1 for exterior penalty = ');
 if c==0
     c = 0.1;              % c to update penalty parameter
 else
     c = 10;
 end
 t = 1;
 x(t,:)= x0;
 R(t)= R0;
 x1 = x(1);
 x2 = x(2);
 x3 = x(3);
 x4 = x(4);
 x5 = x(5);
 x6 = x(6);
 x7 = x(7);
 x8 = x(8);
g1 = (1 - 0.0025 *(x4 + x6));
g2 = (1 - 0.0025*(-x4 +x5 + x7));
g3 = (1-0.01*(-x6 + x8));
g4 = (83333.333 - 833.33252*x4 + x1*x6 - 100*x1);
g5 = (x2*x7 + 1250*x4 - 1250*x5 - x2*x4);
g6 = (x3*x8 - x3*x5 - 2500*x5 - 1250000);
 if (g1 > 0 && g2 >0 && g3 > 0 && g4 >0 && g5 > 0 && g6>0)
     x1 = x(t,1);
       x2 = x(t,2);
       x3 = x(t,3);
       x4 = x(t,4);
       x5 = x(t,5);
       x6 = x(t,6);
       x7 = x(t,7);
       x8 = x(t,8);
       k = 0;
       F = penfun(x(t,:),0);
        P(1,t) = penfun(x(t,:),0);
        G(t,:) = Phase_2(x(t,:),0);
      t = t+1;
      P(1,t) = penfun(G(t-1,:),0);
      err = abs(P(1,t));
     R(1,t) = c*R(1,t-1);
     x(t,:) = G(t-1,:);  
 else 
     for t = 1
        x1 = x(t,1);
        x2 = x(t,2);
        x3 = x(t,3);
        x4 = x(t,4);
        x5 = x(t,5);
        x6 = x(t,6);
        x7 = x(t,7);
        x8 = x(t,8);
          g1 = (1 - 0.0025 *(x4 + x6));
           g2 = (1 - 0.0025*(-x4 +x5 + x7));
           g3 = (1-0.01*(-x6 + x8));
           g4 = (83333.333 - 833.33252*x4 + x1*x6 - 100*x1);
           g5 = (x2*x7 + 1250*x4 - 1250*x5 - x2*x4);
           g6 = (x3*x8 - x3*x5 - 2500*x5 - 1250000);
           F = penfun(x(t,:),0);
          P(1,t) = penfun(x(t,:),R(1,t));
            G(t,:) = Phase_2(x(t,:),R(1,t));     % to perform unconstrained search for given R
             x(t+1,:) = G(t,:);
         t = t+1;
        P(1,t) = penfun(G(t-1,:),R(1,t-1));
         err = abs(P(1,t));
         R(1,t) = c*R(1,t-1);
          x1 = x(t,1);
        x2 = x(t,2);
        x3 = x(t,3);
        x4 = x(t,4);
        x5 = x(t,5);
        x6 = x(t,6);
        x7 = x(t,7);
        x8 = x(t,8);
          g1 = (1 - 0.0025 *(x4 + x6));
            g2 = (1 - 0.0025*(-x4 +x5 + x7));
            g3 = (1-0.01*(-x6 + x8));
            g4 = (83333.333 - 833.33252*x4 + x1*x6 - 100*x1);
            g5 = (x2*x7 + 1250*x4 - 1250*x5 - x2*x4);
            g6 = (x3*x8 - x3*x5 - 2500*x5 - 1250000);
     end
  
 end
 
 a(:,1) = G(t-1,:)';

 while err>tol2
     G(t,:) = Phase_2(x(t,:),R(1,t));
     x(t+1,:) = G(t,:);
     t = t+1;
     R(1,t) = c*R(1,t-1);
     P(1,t) = penfun(G(t-1,:),R(1,t-1));
     err = abs(P(1,t)-P(1,t-1));
       x1 = x(t,1);
        x2 = x(t,2);
        x3 = x(t,3);
        x4 = x(t,4);
        x5 = x(t,5);
        x6 = x(t,6);
        x7 = x(t,7);
        x8 = x(t,8);
          g1 = (1 - 0.0025 *(x4 + x6));
            g2 = (1 - 0.0025*(-x4 +x5 + x7));
            g3 = (1-0.01*(-x6 + x8));
            g4 = (83333.333 - 833.33252*x4 + x1*x6 - 100*x1);
            g5 = (x2*x7 + 1250*x4 - 1250*x5 - x2*x4);
            g6 = (x3*x8 - x3*x5 - 2500*x5 - 1250000);
     a(:,1) = G(t-1,:)';
     if (abs(g1) <tol1 && abs(g2)<tol1 && abs(g3)<tol1 && abs(g4)<tol1 && abs(g5)<tol1 && abs(g6)<tol1)
         break;
     elseif (abs(x(t,1)-x(t-1,1))<10e-5 && abs(x(t,2)-x(t-1,2))<10e-5 && abs(x(t,3)-x(t-1,3))<10e-5 && abs(x(t,4)-x(t-1,4))<10e-5 && abs(x(t,5)-x(t-1,5))<10e-5 && abs(x(t,6)-x(t-1,6))<10e-5 && abs(x(t,7)-x(t-1,7))<10e-5 && abs(x(t,8)-x(t-1,8))<10e-5)
         break;
     end
     
 end
 fval = penfun(x(t,:),0);
fprintf('The minimum point is (%4f,%4f,%4f,%4f,%4f,%4f,%4f,%4f) with a function value of %.4f  \n',a,fval)

%% POWELL CONJUGATE DIRECTION METHOD- PHASE 2
function fin_vec = Phase_2(x0,R)
fprintf('Powell Conjugate Direction Method')
k=1;                                 % iterations' counter
N = length(x0);                      % state's dimension
tol3 = 1e-4;                          % tolerance to stop the algorithm 
condition = true;
while condition
    s1 = eye(N,N);
    x = zeros(N,N);
    for j = 1:N
        alpha = PhaseI(x0,s1,tol3,j,R);
        for i =1:N
            x(j,:) = x0 + alpha*s1(j,:);

        end
        x0 = x(j,:);
    end
    % NEW POINT
    alpha = PhaseI(x0,s1,tol3,j,R);
    x1 = x0 + alpha*s1(1,:);
    X = x1;
    Y = x(1,:);
     % UPDATE and STORAGE
    d = X-Y;
    p = d'*d;
    q = sqrt(p);
    norm_vec = d/q;

    fin_vec = X;
    if (q<tol3)
        break;
    elseif(k>300)
        break;
    else
        s1(1,:) = norm_vec;
        k = k+1;
    end
end
fprintf('\n*************************\n');
fprintf('\nTotal number of Iterations for POWELL METHOD: %d\n', k);
fprintf("The minimum function value is :%f \n", penfun(x,R));
fmt=['The minimum solutions is =' repmat(' %1.0f',1,numel(X))];
fprintf(fmt,X)
end
%% Objective function
function P = penfun(x,R)
     x1 = x(1);
     x2 = x(2);
     x3 = x(3);
     x4 = x(4);
     x5 = x(5);
     x6 = x(6);
     x7 = x(7);
     x8 = x(8);
     P = (x1 + x2 +x3);
     g1 = (1 - 0.0025 *(x4 + x6));
     g2 = (1 - 0.0025*(-x4 +x5 + x7));
     g3 = (1-0.01*(-x6 + x8));
     g4 = (83333.333 - 833.33252*x4 + x1*x6 - 100*x1);
     g5 = (x2*x7 + 1250*x4 - 1250*x5 - x2*x4);
     g6 = (x3*x8 - x3*x5 - 2500*x5 - 1250000);
     if g1<0
         P = P + R*g1^2;
     end
     if g2 <0
         P = P + R*g2^2;
     end
     if g3 < 0
         P = P + R*g3^2;
     end
     if g4 < 0
         P = P + R*g4^2;
     end
     if g5 < 0
         P = P + R*g5^2;
     end
     if g6 < 0
         P = P + R*g6^2;
     end
     
end
%% MAPPING FUNCTION (FOR UNIDIRECTIONAL SEARCH)
function fun_val = Obj_Fun(y,x0,s1,j,R)

x0 = x0 + y*s1(j,:);
fun_val = penfun(x0,R);
end

%% First Order Derivative Function
function fun_diff=Diff_fun(y, x0,s1,j,R)
h=0.001;%%Step size
fun_diff=(Obj_Fun(y+h,x0,s1,j,R)-Obj_Fun(y-h,x0,s1,j,R))/(2*h);
end
%% MAIN PROGRAM FOR PHASE 1
function z=PhaseI(x0,s1,tol,j,R)
[p,q]=bounding_phase(x0,s1,j,R);
[z] = secant_method(x0,s1, p, q, tol,j,R);
end

%% BOUNDING PHASE METHOD - BRACKETING METHOD
function [p,q]=bounding_phase(x0,s1,j,R)
fprintf('Bounding Phase Method')
delta = 0.04;
y0 = 100;
% y0 =2.6;
feval = 0; % Number of function evaluations
y1=y0-delta;y2=y0+delta;
f1 = Obj_Fun(y1,x0,s1,j,R);
f2 = Obj_Fun(y0,x0,s1,j,R);
f3 = Obj_Fun(y2,x0,s1,j,R);
feval = feval + 3;
out = fopen('Bounding_Phase_Method.out', 'w'); % Output file
fprintf(out, '#It\tx1\tx2\tx3\tf(x1)\tf(x2)\tf(x3)\n');
condition=true;
if (f1<f2 && f2<f3) %% Minimization
    delta=-abs(delta);

elseif f3<=f2 && f2<=f1 %% Minimization
    delta= abs(delta);
    y1 = y2;
    f1 = f3;
end
i=0;yd=y0;
while condition
    fprintf(out, '%d\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\n',i,y1,y0,y2,f1,f2,f3);
    q=y0+(2^i)*delta;
    fq=Obj_Fun(q,x0,s1,j,R);
    feval = feval + 1;

    if fq<f2 %% Minimization
        yd=y0;
        y0=q;
        f2=fq;
        i=i+1;
        
    else
        p=yd;
        condition=false;
    end
    
end
fprintf('\n*************************\n');
    fprintf('The minimum point lies between (%8.3f, %8.3f)', y1, y2);
    fprintf('\nTotal number of function evaluations: %d\n', feval)
 % Store in the file
    fprintf(out, '\nThe minimum point lies between (%8.3f, %8.3f)', y1, y2);
    fprintf(out, '\nTotal number of function evaluations: %d', feval);

fclose(out);
end
%% SECANT METHOD - GRADIENT BASED METHOD
function [z] = secant_method(x0,s1, p, q, tol,j,R)
fprintf('Secant Method')

x1=p; x2=q; feval =0;
%1st derivative
fdx1=Diff_fun(x1,x0,s1,j,R);
fdx2=Diff_fun(x2,x0,s1,j,R);
feval = feval + 4;
out = fopen('Secant Method.out', 'w'); % Output file
fprintf(out, '#It\tx1\tx2\tx3\tf(x1)\tf(x2)\tf(x3)\n');
i=0;fdz=1;
while abs(fdz)>tol
    
    z=x2-fdx2*(x2-x1)/(fdx2-fdx1);
    fdz=Diff_fun(z,x0,s1,j,R);
    feval = feval + 2;
    fprintf(out, '%d\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\n',i,x1,x2,z,fdx1,fdx2,fdz);
    if fdz<0
        x1=z;
        fdx1=fdz;
    elseif fdz>0
        x2=z;
        fdx2=fdz;
    end
    
    i=i+1;
    
end
fprintf('\n*************************\n');
    fprintf('The minimum point is %d', z);
    fprintf('\nTotal number of function evaluations: %d\n', feval);
    
    % Store in the file
    fprintf(out, '\nThe minimum point is %d', z);
    fprintf(out, '\nTotal number of function evaluations: %d', feval);
fclose(out);
end