x=input('Is equation is a polynomial (Y/N) \n','s');
clc
if ~strcmp(x,'Y')
    fprintf('(1) Bisection \n(2) False-position \n(3) Fixed-Point \n(4) Newton-Raphson \n(5) Secant \n');
    y=input('');
    clc
    if y==1
        i=0;
        flag=0;
        f_x = input('Give an equation in x:\n','s');
        d=input('Enter two starting points using [a b]\n');
        fprintf('Stopping criteria:\n(i) Convergence criterion for relative approximate errors in successive iterations. The approximate relative error may be for the root for some methods, interval for some methods and factor values for some methods. \n(ii) Convergence criteria for the function value, i.e., how close f(x) is to zero. \n(iii) Maximum iteration number (an integer)\n');
        fprintf('Enter maximum relative approximate error and number of iterations line by line\n');
        tol=input('');
        n=input('');
        b=d(2);
        a=d(1);
        x1=[];
        y1=[];
        %disp(feval(inline(f_x,'x'),a));
        %disp(feval(inline(f_x,'x'),b));
        while abs(b-a)>=tol/100
            c=(a+b)/2;
            i=i+1;
            x1=[x1,i];
            y1=[y1,abs(b-a)];
            fc=feval(inline(f_x,'x'),c);
            %disp(fc);
            %disp(abs(b-a));
            if feval(inline(f_x,'x'),a)*fc<0
                b=c;
                fb=fc;
            elseif feval(inline(f_x,'x'),b)*fc<0
                a=c;
                fa=fc;
            else
                flag=2;
                break
            end
            if i==n
                flag=3;
                break;
            end
            if abs(b-a)>tol
                flag=1;
            end
        end
        disp(c);
        disp(flag);
        x = -20:pi/100:20;
        y = feval(inline(f_x,'x'),x);
        plot(x,y);
        figure;plot(x1,y1);
    end
    if y==2
        i=0;
        flag=0;
        f_x = input('Give an equation in x:\n','s');
        d=input('Enter two starting points using [a b]\n');
        fprintf('Stopping criteria:\n(i) Convergence criterion for relative approximate errors in successive iterations. The approximate relative error may be for the root for some methods, interval for some methods and factor values for some methods. \n(ii) Convergence criteria for the function value, i.e., how close f(x) is to zero. \n(iii) Maximum iteration number (an integer)\n');
        fprintf('Enter maximum relative approximate error and number of iterations line by line\n');
        tol=input('');
        n=input('');
        b=d(2);
        a=d(1);
        x1=[];
        y1=[];
        while abs(b-a)>=tol/100
            i=i+1;
            c=a-((b-a)/(feval(inline(f_x,'x'),b)-feval(inline(f_x,'x'),a)))*feval(inline(f_x,'x'),a);
            x1=[x1,i];
            y1=[y1,abs(b-a)];
            fc=feval(inline(f_x,'x'),c);
            %disp(feval(inline(f_x,'x'),a)*fc);
            if feval(inline(f_x,'x'),a)*fc<0
                b=c;
                fb=fc;
            elseif feval(inline(f_x,'x'),b)*fc<0
                a=c;
                fa=fc;
            else
                break
            end
            if i==n
                flag=3;
                break
            end
            if abs(b-a)>tol
                flag=1;
            end
        end
        disp(c);
        disp(flag);
        x = -20:pi/100:20;
        y = feval(inline(f_x,'x'),x);
        plot(x,y);
        figure;plot(x1,y1);
    end
    if y==3
        i=0;
        flag=0;
        f_x = input('Give an equation in x:\n','s');
        g_x=input('Enter g(x) such that f(x) can be rearranged as x=g(x) and one starting point.\n','s');
        d=input('');
        fprintf('Stopping criteria:\n(i) Convergence criterion for relative approximate errors in successive iterations. The approximate relative error may be for the root for some methods, interval for some methods and factor values for some methods. \n(ii) Convergence criteria for the function value, i.e., how close f(x) is to zero. \n(iii) Maximum iteration number (an integer)\n');
        fprintf('Enter maximum relative approximate error and number of iterations line by line\n');
        tol=input('');
        n=input('');
        a=d;
        x1=[];
        y1=[];
        b=feval(inline(g_x,'x'),a);
        %disp(feval(inline(f_x,'x'),a));
        %disp(b);
        while abs((b-a)/a)>tol/100
            i=i+1;
            a=b;
            b=feval(inline(g_x,'x'),a);
            i=i+1;
            fc=feval(inline(f_x,'x'),b);
            x1=[x1,i];
            y1=[y1,abs((b-a)/a)];
            %disp(fc);
            %disp(abs((b-a)/b));
            if feval(inline(f_x,'x'),a)==0
                flag=2;
                break
            end
            if i==n
                flag=3;
                break
            end
            if abs((b-a)/a)>tol
                flag=1;
                
            end
        end
        disp(a);
        disp(flag);
        x = -20:pi/100:20;
        y = feval(inline(f_x,'x'),x);
        plot(x,y);
        figure;plot(x1,y1);
    end
    if y==4
        i=0;
        flag=0;
        f_x = input('Give an equation in x:\n','s');
        g_x=input('Enter f''(x) and one starting point.\n','s');
        d=input('');
        fprintf('Stopping criteria:\n(i) Convergence criterion for relative approximate errors in successive iterations. The approximate relative error may be for the root for some methods, interval for some methods and factor values for some methods. \n(ii) Convergence criteria for the function value, i.e., how close f(x) is to zero. \n(iii) Maximum iteration number (an integer)\n');
        fprintf('Enter maximum relative approximate error and number of iterations line by line\n');
        tol=input('');
        n=input('');
        a=d;
        x1=[];
        y1=[];
        b=a-(feval(inline(f_x,'x'),a)/feval(inline(g_x,'x'),a));
        %disp(feval(inline(f_x,'x'),a));
        %disp(b);
        while abs((b-a)/a)>tol/100
            i=i+1;
            a=b;
            b=a-(feval(inline(f_x,'x'),a)/feval(inline(g_x,'x'),a));
            i=i+1;
            fc=feval(inline(f_x,'x'),b);
            x1=[x1,i];
            y1=[y1,abs((b-a)/a)];
            %disp(fc);
            %disp(abs((b-a)/b));
            if feval(inline(f_x,'x'),a)==0
                flag=2;
                break
            end
            if i==n
                flag=3;
                break
            end
            if abs((b-a)/b)>tol
                flag=1;
                
            end
        end
        disp(a);
        disp(flag);
        x = -20:pi/100:20;
        y = feval(inline(f_x,'x'),x);
        plot(x,y);
        figure;plot(x1,y1);
    end
    if y==5
        i=0;
        flag=0;
        f_x = input('Give an equation in x:\n','s');
        d=input('Enter two starting points using [a b]\n');
        fprintf('Stopping criteria:\n(i) Convergence criterion for relative approximate errors in successive iterations. The approximate relative error may be for the root for some methods, interval for some methods and factor values for some methods. \n(ii) Convergence criteria for the function value, i.e., how close f(x) is to zero. \n(iii) Maximum iteration number (an integer)\n');
        fprintf('Enter maximum relative approximate error and number of iterations line by line\n');
        tol=input('');
        n=input('');
        x1=[];
        y1=[];
        x_prev=d(1);
        x_current=d(2);
        x_next=x_current-((x_current-x_prev)/(feval(inline(f_x,'x'),x_current)-feval(inline(f_x,'x'),x_prev)))*feval(inline(f_x,'x'),x_current);
        while abs((x_next-x_current)/x_current)>=tol/100
            i=i+1;
             x1=[x1,i];
            y1=[y1,abs((x_next-x_current)/x_current)];
            x_prev=x_current;
            x_current=x_next;
            x_next=x_current-((x_current-x_prev)/(feval(inline(f_x,'x'),x_current)-feval(inline(f_x,'x'),x_prev)))*feval(inline(f_x,'x'),x_current);
            fc=feval(inline(f_x,'x'),x_next);
            %disp(feval(inline(f_x,'x'),a)*fc);
            if feval(inline(f_x,'x'),a)*fc<0
                b=c;
                fb=fc;
            elseif feval(inline(f_x,'x'),b)*fc<0
                a=c;
                fa=fc;
            else
                break
            end
            if i==n
                flag=3;
                break
            end
            if abs((x_next-x_current)/x_current)>=tol
                flag=1;
            end
        end
        disp(x_next);
        disp(flag);
        x = -20:pi/100:20;
        y = feval(inline(f_x,'x'),x);
        plot(x,y);
        figure;plot(x1,y1);
    end
else
    s=input('Enter the function ','s');  %% taking input from the user
    s1='@(x)';
    c=strcat(s1,s);
    f=str2func(c);
    %n=input('Enter the degree of the polynomial\n');
    %a=input('Enter the coefficients as [a0 a1 a2 .... an]\n');
    fprintf('(1) Muller \n(2) Bairstow \n');
    y=input('');
    if y==1
        disp('Enter starting value x1');
        x_ppprev=input('');
        disp('Enter starting value x2');
        x_pprev=input('');
        disp('Enter starting value x3');
        x_prev=input('');

        disp('Enter max no of iterations ');
        max_no_interation=input('');
        disp('Enter max_relative_error ');
        max_relative_error=input('');
        x1=[];
        y1=[];
        w=2;
        rel_error=100*ones(max_no_interation,1);
        flag=0;
        while(abs(rel_error(i-1))>max_relative_error & i-2< max_no_interation)
            a=((f(x_prev)-f(x_pprev))/(x_prev-x_pprev)- (f(x_pprev)-f(x_ppprev))/(x_pprev-x_ppprev))/(x_prev-x_ppprev);
            b=(f(x_prev)-f(x_pprev))/(x_prev-x_pprev)+a*(x_prev-x_pprev);
            c=f(x_prev);
            i=i+1;
            x1=[x1,i];
            y1=[y1,abs(rel_error(i-1))];
            if(b<0)
            delta_x=-2*c/(b-sqrt(b^2-4*a*c));
            else 
            delta_x=-2*c/(b+sqrt(b^2-4*a*c));
            end
    
            x_next=x_prev+delta_x;
            if(i-3>0)
            rel_error(i-1)=((x_next-x_prev)*100)/x_next;
            end 
            if(abs(rel_error(i-1))>max_relative_error)
                flag=1;
            end
            if(w==max_no_interation)
                flag=3;
            end
            x_ppprev=x_pprev;
            x_pprev=x_prev;
            x_prev=x_next;
        end
       plot(x1,y1);
       disp(x_prev);
       disp(flag);
    end
    if y==2;
        
    coff=input('Enter the coeff in Increasing Order in []:');
    r=input('Enter s :');
    s=input('Enter r :');
    maxit=input('Enter Maximum no.of Iterations :');
    errmax=input('Enter the Maximum Relative approximate error in(%) :');
    degree=length(coff)-1;
	scoff=coff;
	sdegree=degree;
    arrc=[];
    d=degree;
    arrb=[];
    flag=0;
    count=0;
    yr=[];
    ys=[];
    iee=0;
    x1=[];
    while(degree>0)
    iee=iee+1;
    x1=[x1,iee];
    for j=1:maxit
	
    d=degree;
    arrb(d+1)=coff(d+1);
    arrb(d)=coff(d)+r*arrb(d+1);
    d=d-2;
    
    while(d>=0)
          arrb(d+1)=coff(d+1) + r*arrb(d+2) + s*arrb(d+3);
          d=d-1; 
    end
    
    d=degree;
    arrc(d+1)=arrb(d+1);
    arrc(d)=arrb(d)+ r*arrc(d+1);    
    d=d-2;
    
     while(d>=0)
          arrc(d+1)=arrb(d+1) + r*arrc(d+2) + s*arrc(d+3);
          d=d-1; 
     end
     
     ds = ( arrb(1)*arrc(3) - arrb(2)*arrc(2) )/( arrc(4)*arrc(2) - arrc(3)*arrc(3) );
     dr = ( arrb(1)*arrc(4) - arrb(2)*arrc(3) )/( arrc(3)*arrc(3) - arrc(2)*arrc(4) );
     r=r+dr;
     s=s+ds;
     
     err_r = abs(dr/r)*100 ;
     err_s = abs(ds/s)*100 ;
     yr=[yr,abs(dr/r)*100];
     ys=[ys,abs(ds/s)*100];
     if ( errmax > err_s ||  errmax > err_r )
         
         root1 = ( r + sqrt(r*r + 4*s) )/2;
         root2 = ( r - sqrt(r*r + 4*s) )/2;
         fprintf('Solution is ');
         
         if(count==0)
             fprintf(' %f %f ',root1,root2);
         end
         count=1;
         flag=1;
         break; 
     end
	 
    end
	
    if flag==0
        break;
    end
	
     for i=1:degree-1
         coff(i)=arrb(i+2);
     end
	 
     degree=degree-2;
	 
     if degree == 2
         root1 = ( -(arrb(4)) + sqrt(arrb(4)*arrb(4) - 4*arrb(3)*arrb(5)) )/(2*arrb(5));
         root2 = ( -(arrb(4)) - sqrt( arrb(4)*arrb(4) - 4*arrb(3)*arrb(5) ))/(2*arrb(5));
         fprintf('%f %f\n ',root1,root2);
         break;
		 
     elseif degree == 1
         root = -arrb(3)/arrb(4);
         fprintf('%f\n ',root);
         break;
		 
     end      
         
    end 
	
	f=zeros(1,101);
   
	
    for k=-50:50
	
        for i=1:sdegree+1
	
            f(k+51)=f(k+51)+scoff(i)*((k)^(i-1)); %x^4-7.4*x^3+20.44*x^2-24.184*x+9.6448  [9.6448 -24.184 20.44 -7.4 1]
	 
        end
	   
    end
	   
	 disp(iee);
 figure;plot(x1,yr);
   
    
    end   
end
