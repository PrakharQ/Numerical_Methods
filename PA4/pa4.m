disp('Enter Grid Size');
h=input('');
p=@(x)1;
q=@(x)(-(x+3)/(x+1));
r=@(x)((x+3)/(x+1)^2);
s=@(x)(2*(x+1)+3*r(x));
l=[];
d=[];
u=[];
b=[];
x=0;
z=1;
for i=1:2/h-1
    x=x+h;
    l(z)=(p(x)/(h*h))-(q(x)/(2*h));
    d(z)=(-2*p(x)/(h*h))+r(x);
    u(z)=(p(x)/(h*h))+(q(x)/(2*h));
    b(z)=s(x);
    z=z+1;
end
%disp(l);
ch=input('1)2nd order Backward Difference\n2)2nd order Central Difference with Ghost Node\n','s');
if(ch=='1')
    b(1)=b(1)-l(1)*5;
    l(z-1)=l(z-1)-u(z-1)/3;
    d(z-1)=d(z-1)+4*u(z-1)/3;
    %l(1)=0;
    %u(z)=0;
    n = length(b);
    v = zeros(n,1);   
    y = v;
    w = d(1);
    y(1) = b(1)/w;
    for i=2:n
    v(i-1) = u(i-1)/w;
    w = d(i) - l(i)*v(i-1);
    y(i) = ( b(i) - l(i)*y(i-1) )/w;
    end
    for j=n-1:-1:1
    y(j) = y(j) - v(j)*y(j+1);
    end
    disp(y);
    x=h:h:2-h;
    plot(x,y);
end
if(ch=='2')
    x=x+h;
    l(z)=(p(x)/(h*h))-(q(x)/(2*h));
    d(z)=(-2*p(x)/(h*h))+r(x);
    u(z)=(p(x)/(h*h))+(q(x)/(2*h));
    b(z)=s(x);
    b(1)=b(1)-l(1)*5;
    l(z)=l(z)+u(z);
    n = length(b);
    v = zeros(n,1);   
    y = v;
    w = d(1);
    y(1) = b(1)/w;
    for i=2:n
    v(i-1) = u(i-1)/w;
    w = d(i) - l(i)*v(i-1);
    y(i) = ( b(i) - l(i)*y(i-1) )/w;
    end
    for j=n-1:-1:1
    y(j) = y(j) - v(j)*y(j+1);
    end
    disp(y);
    x=h:h:2;
    plot(x,y);
end
