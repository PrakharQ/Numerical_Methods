%% input from file
filename=input('Enter the name of the file from where to take inputs\n','s');
fileid=fopen(filename,'r');  
formatSpec = '%f';
format long;
n=fscanf(fileid,formatSpec,1);
x=fscanf(fileid,formatSpec,[1 n]); 
y=fscanf(fileid,formatSpec,[1 n]); 
%%
disp('For Least Square Polynomial Choose A ');
disp('For Langrage Interpolation Choose B ');
disp('For use of Cubic Splines Choose C');
choice = input('What you want to do ?\n', 's');
if choice == 'A'
m = input('Specify the Order of Approximate Polynomial\n');
choice = input('Specify Whether intercept of Polynomials is Zero[Y] or not[N]\n','s');
if choice == 'N'
    N = zeros(m+1);
    c = zeros(m+1,1);
    
    for i = 1:m+1
        for j = i:m+1
            for k =1:n
                N(i,j) = N(i,j) + power(x(k),i+j-2) ;
            end
            N(j,i) = N(i,j);
        end
    end
    for j =1:m+1
        for k =1:n
            c(j) = c(j) + power(x(k),j-1)*y(k); 
        end
    end
    coeff = N\c ;
     y0 = [];
    for x0 = x(1):0.1:x(n)
        temp = 0;
        for j = 1:m+1
         temp = temp + power(x0,j-1)*coeff(j) ;
        end
        y0 = [y0,temp] ;
    end
    x0 = x(1):0.1:x(n) ;
    plot(x0,y0,'b',x,y,'or') ;
    grid on;
    y_cal =  zeros(1,n);
    for i=1:n
        for j=1:m+1
            y_cal(i) = y_cal(i) +  power(x(i),j-1)*coeff(j) ;
        end
    end
    r_2 = (sum((y-mean(y_cal)).^2)-sum((y-y_cal).^2))/sum((y-mean(y_cal)).^2) ;
    r_2=abs(r_2);
    
     %% Output to File
     filename = input('Enter The Name of Output file\n','s') ;
     fileID = fopen(filename,'w') ;
     
     for i = 0:1:m
         fprintf(fileID,'The Coefficients are a%d \t %4.4f \n',i,coeff(i+1)) ;
     end
     
        fprintf(fileID,'The value or r_square is \t %f',r_2) ;
     fclose(fileID) ;
else
    N = zeros(m);
    c = zeros(m,1);
    for i = 1:m
        for j = i:m
            for k =1:n
                N(i,j) = N(i,j) + power(x(k),i+j) ;
            end
            
            N(j,i) = N(i,j);
        end
    end
    for j =1:m
        for k =1:n
            c(j) = c(j) + power(x(k),j)*y(k); 
        end 
    end
    coeff = N\c ;
    y0 = [];
    for x0 = x(1):0.1:x(n)
        temp = 0;
        for j = 1:m
         temp = temp + power(x0,j)*coeff(j) ;
        end
        y0 = [y0,temp] ;
    end
    x0 = x(1):0.1:x(n) ;
    plot(x0,y0,'b',x,y,'or') ;
    grid on;
    y_cal =  zeros(1,n);
    for i=1:n
        for j=1:m
            y_cal(i) = y_cal(i) +  power(x(i),j)*coeff(j) ;
        end
    end
    r_2 = (sum((y-mean(y_cal)).^2)-sum((y-y_cal).^2))/sum((y-mean(y_cal)).^2) ;
    r_2=abs(r_2);
    
    %% Output to File
    
     filename = input('Enter The Name of Output file\n','s') ;
     fileID = fopen(filename,'w') ;
     
     for i = 1:m
         fprintf(fileID,'The Coefficients are a%d \t %4.4f \n',i,coeff(i)) ;
     end
     
        fprintf(fileID,'The value or r_square is \t %f',r_2) ;
     fclose(fileID) ;
    %%
    
end
elseif choice == 'B'
     syms x_var;
f=@(x_var)0;
Li=@(x_var)1;
for i=1:n
    for j=1:n
        if(i~=j)
            Li=@(x_var)Li(x_var)*(x_var - x(j))/(x(i)-x(j));
            
        end
    end
    f=@(x_var)f(x_var) + Li(x_var)*y(i);
    Li=@(x_var)1;
end
l=1;
for p=x(1):0.01:x(n)
    w(l)=p;
    a=f(p);
    z(l)=a;
    l=l+1;
end
scatter(x,y);
hold on
plot(w,z);
    
else
   dv_df(n-1) = 0;
   h(n-1) = 0 ;
   col = zeros(n-2,1) ;
   for i =1:n-1
       dv_df(i) = (y(i+1)-y(i))/(x(i+1)-x(i)) ;
       h(i) = (x(i+1)-x(i)) ;
   end
   
   H = zeros(n-2,n) ;
   A = zeros(n-2) ;
   sol = [] ;
   for i =1:n-2
       H(i,i) = h(i) ;
       H(i,i+1) = (h(i)+h(i+1))*2 ;
       H(i,i+2) = h(i+1) ;
       col(i) = 6*(dv_df(i+1)-dv_df(i)) ;
   end
   disp('from the Given Option Choose');
   disp('1. Natural Spline');
   disp('2. Not-a-Knot');
   disp('3. Periodic');
   disp('4. Clamped Spline');
   
   choice = input('Enter Your Choice\n','s');
   switch(choice)
       case '1'
                A(:,1) =  H(:,2) ;
                A(:,n-2) = H(:,n-1) ;
                A(:,2:n-3) = H(:,3:n-2) ;
                v0 = 0 ; vn = 0;
                
                sol = A\col ;
                sol = [v0;sol;vn] ;
       case '2'
                A(:,1) = ((h(2)+h(1))/h(2))*H(:,1) + H(:,2) ;
                A(:,2) = -1*(h(1)/h(2))*H(:,1) + H(:,3) ;
                A(:,n-2) = H(:,n-1) + ((h(n-1)+h(n-2))/h(n-2))*H(:,n) ;
                A(:,n-3) = H(:,n-2) -(h(n-1)/h(n-2))*H(:,n) ;
                A(:,3:n-4) = H(:,4:n-3) ;
                
                sol = A\col ;
                v0 = sol(1) - (sol(2)-sol(1))*(h(1)/h(2)) ;
                vn = sol(n-2) + (sol(n-2)-sol(n-3))*(h(n-1)/h(n-2)) ;
                
                sol = [v0;sol;vn] ;
       case '3'
                A(:,1) = H(:,2) + H(:,n) ;
                A(:,n-2) = H(:,n-1) + H(:,1) ;
                A(:,2:n-3) = H(:,3:n-2) ;
                
                sol = A\col ;
                v0 = sol(n-2);
                vn = sol(1);
                
                sol = [v0;sol;vn] ;
       case '4'
                u0 = input('Enter Derivative at First Node \n') ;
                un = input('Enter Derivative at Last Node \n') ;
                col2 = col - ((3*(dv_df(1)-u0))/h(1))*H(:,1) ;
                A(:,1) = H(:,2) - (H(:,1)/2) ;
                A(:,n-2) = H(:,n-1) - (H(:,n)/2) ;
                A(:,2:n-3) = H(:,3:n-2) ;
                col2 = col2 -((3*(un - dv_df(n-1)))/h(n-1))*H(:,n) ;
                
                sol = A\col2 ;
                
                v0 = (((6*(dv_df(1)-u0))/h(1))-sol(1))/2 ;
                vn = (((6*(un - dv_df(n-1)))/h(n-1))-sol(n-2))/2 ;
                
                sol = [v0;sol;vn] ;
   end
   for i = 1:n-1
       y_val = zeros(1,(x(i+1)-x(i))*100 + 1) ;
       k = 1 ;
       for x_val = x(i):0.01:x(i+1)
       f = (sol(i+1)/6)*(( power(x_val-x(i),3))/h(i) - h(i)*(x_val-x(i))) + (sol(i)/6)*( (-power(x_val-x(i+1),3))/h(i) + h(i)*(x_val-x(i+1))) + (y(i+1)/h(i))*(x_val-x(i)) - (y(i)/h(i))*(x_val-x(i+1)) ;
       y_val(k) = f  ;
       k = k+1 ;
       end
   x_val = x(i):0.01:x(i+1) ;
   plot(x_val,y_val,'b') ;    
   hold on
   end
   grid on
   plot(x,y,'or')
   hold off
   %% Output to File
    
     filename = input('Enter The Name of Output file\n','s') ;
     fileID = fopen(filename,'w') ;
     
     for i = 1:n-1
         a = (sol(i+1)-sol(i))/(6*h(i)) ;
         b = (sol(i+1)*x(i)-sol(i)*x(i+1)) /(2*h(i)) ;
         c_ = sol(i)*(h(i)/6 - x(i+1)^2/(2*h(i))) + sol(i+1)*(-h(i)/6 + x(i)^2/(2*h(i))) + dv_df(i) ;
         d =  sol(i+1)*(-x(i)^3/h(i) + h(i)*x(i))/6 + sol(i)*(x(i+1)^3/h(i) - h(i)*x(i+1))/6 + (y(i)*x(i+1)-y(i+1)*x(i))/h(i) ;
         
         fprintf(fileID,'The Interval for Interpolation %f to %f \n',x(i),x(i+1)) ;
         fprintf(fileID,'The Coefficients are a3 = %4.4f \t a2 = %4.4f \t a1 = %4.4f \t a0 = %4.4f \n',a,b,c_,d) ;
         fprintf(fileID,'The Value of the first derivative at  first node is %4.4f \t and at second node is %4.4f \n',3*a*x(i)^2 + 2*b*x(i) + c_,3*a*x(i+1)^2 + 2*b*x(i+1) + c_)  ;
         fprintf(fileID,'The Value of the second derivative at  first node is %4.4f \t and at second node is %4.4f \n',sol(i),sol(i+1)) ;
     
     end
         fclose(fileID) ;    
         
end