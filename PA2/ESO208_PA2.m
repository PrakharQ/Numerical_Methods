x=input('A. Solve a System of Equation \nB. Perform a LU decomposition \nC. Perform a Matrix Inversion\n','s');
clc
if strcmp(x,'A')
    y=input('Is the system tri-diagonal (Y/N)\n','s');
    if strcmp(y,'Y')
        filename=input('Enter the name of the file from where to take inputs\n','s');
        fileID = fopen(filename,'r');
        formatSpec = '%f';
        N=fscanf(fileID,formatSpec,1);
        A = fscanf(fileID,formatSpec,[4,4]);
        A=transpose(A);
        a=A(1:1,:);
        b=A(2:2,:);
        c=A(3:3,:);
        r=A(4:4,:);
        
        %[a,b,c,r,N] = textscan(filename,'%f %f %f %f %f',1)
        beta = b(1) ;
        u(1) = r(1)/beta ;
        for j = 2:N
        gamma(j) = c(j-1)/beta ;
        beta = b(j)-a(j)*gamma(j) ;
        if (beta==0)
        fprintf(1,'The tridiagonal solver failed...'); 
        pause
        end
        u(j) = (r(j)-a(j)*u(j-1))/beta ;
        end
        for j = 1:(N-1)
        k = N-j ;
        u(k) = u(k) - gamma(k+1)*u(k+1) ;
        end
        %fileID = fopen('output1.txt','s');
        %fprintf(fileID,'%f',A);
        fileID = fopen('outputQ1_Thomas.txt','w');
        fprintf(fileID,'Solution=\n%f  \n  \n',u);
        fclose(fileID);
        type outputQ1_Thomas.txt
        disp(u);
    else
        filename=input('Enter the name of the file from where to take inputs\n','s');
        fileID = fopen(filename,'r');
        formatSpec = '%f';
        A = fscanf(fileID,formatSpec,[3,4]);
        C=A;
        n=3;
        %B = fscanf(fileID,formatSpec,[3,4]);
        %B=transpose(B);
        k=1;
        maxi=n;
        while(k<n)
  
        for i=k:n-1      %find max in column
        if(abs(C(i,k))<=abs(C(i+1,k)))
            maxi=i+1;
        end
        end
        t=C(k,:);%bring  pivot row up
        C(k,:)=C(maxi,:);
        C(maxi,:)=t;
        for l=k+1:n
        if(C(k,k)~=0)
            C(l,:)=C(l,:)-C(l,k)*C(k,:)/C(k,k);
        end
        end
        k=k+1;
        end
        x=(inv(C(1:3,1:3))*C(:,4));
        fileID = fopen('outputQ1_GEpartial.txt','w');
        fprintf(fileID,'Solution=\n%f  \n  \n',x);
        fclose(fileID);
        type outputQ1_GEPartial.txt
        disp(x);
    end
end  
if strcmp(x,'B')
    y=input('Is the matrix symmetric and positive definite>(Y/N)','s');
    if strcmp(y,'Y')
        filename=input('Enter the name of the file from where to take inputs\n','s');
        fileID = fopen(filename,'r');
        formatSpec = '%f';
        n= fscanf(fileID,formatSpec,1);
        A = fscanf(fileID,formatSpec,[n,n]);
        L = zeros(n);
        for i=1:n
        L(i,i) = sqrt(A(i,i) - L(i,:)*L(i,:)');
        for j=(i + 1):n
        L(j,i) = (A(j,i) - L(i,:)*L(j,:)')/L(i,i);
        end
        end
        disp(L);
        U=L';
        disp(U);
        fileID = fopen('outputQ1_LU1.txt','w');
        fprintf(fileID,'%f    ',L(3:1));
        fprintf(fileID,'\n ');
        fprintf(fileID,'%f    ',L(3:1));
        fprintf(fileID,'%f    ',L(3:1));
        fprintf(fileID,'%f    ',U);
        fclose(fileID);
        type outputQ1_LU1.txt
    end
    if strcmp(y,'N')
        filename=input('Enter the name of the file from where to take inputs\n','s');
        fileID = fopen(filename,'r');
        formatSpec = '%f';
        n= fscanf(fileID,formatSpec,1);
        A = fscanf(fileID,formatSpec,[n,n]);
        w=input('Doolittle(D) or Crout(C) decomposition\n','s');
        if strcmp(w,'D')
            L=zeros(n,n);
            U=zeros(n,n);
            for k=1:n 
            U(1,k)=A(1,k);
            end
            for k=1:n 
            L(k,k)=1;
            end
 
            k=1;
            for i=2:n
            for k=1:i-1
            temp1=A(i,k);
                if(k==1 )
                    p=1;
                    if( U(p,k)~=0)
                        L(i,p)=A(i,k)/U(p,k);
                    end
                else
                     for p=1:k-1
                        temp1=temp1-L(i,p)*U(p,k);
                     end
                end
            
            L(i,k)=temp1/U(k,k);
            end;
            for k=i:n
               temp=A(i,k);
                for p=1:i-1
                    temp=temp-L(i,p)*U(p,k);
                end
            U(i,k)=temp;
            end
            end
            disp(L);
            disp(U);
        end
        
        if strcmp(w,'C')
            L=zeros(n,n);
            U=zeros(n,n);
            for k=1:n 
            L(k,1)=A(k,1);
            end
            for k=1:n 
            U(k,k)=1;
            end
 
            for j=1:n
            if( L(1,1)~=0)
            U(1,j)=A(1,j)/L(1,1);
            end
            end
 
            for i=2:n
            for j=2:i
            temp=A(i,j);
                for p=1:j-1
                    temp=temp-L(i,p)*U(p,j);
                end
              L(i,j)=temp;
            end
            for j=i:n
              temp1=A(i,j);
                   for p=1:i-1
                        temp1=temp1-L(i,p)*U(p,j);
                    end
                      U(i,j)=temp1/L(i,i);
            end       
            end;
            disp(L);
            disp(U);
            fileID = fopen('outputQ1_LU2.txt','w');
        fprintf(fileID,'Solution=\n%f  \n  \n',L);
        fprintf(fileID,'Solution=\n%f  \n  \n',U);
        fclose(fileID);
        type outputQ1_LU2.txt
        end
   end
        
end
if strcmp(x,'C')
    filename=input('Enter the name of the file from where to take inputs\n','s');
        fileID = fopen(filename,'r');
        formatSpec = '%f';
        A = fscanf(fileID,formatSpec,[4,4]);
        [m,n] = size(A);
        x = 1; % iterator for elimination matrix 3rd dimension
        for j = 1:1:n-1,
    for i = j+1:1:m,
        E(:,:,x) = eye(n);
        E(i,j,x) = -A(i,j) / A(j,j);
        A = E(:,:,x) * A;
        x = x + 1;
    end
        end
        x = 1;
        for j = n:-1:2,
    for i = j-1:-1:1,
        P(:,:,x) = eye(n);
        P(i,j,x) = -A(i,j) / A(j,j);
        A = P(:,:,x) * A;
        x = x + 1;
    end
        end
        for i = 1:1:n,
    S(:,:,i) = eye(n);
    S(i,i,i) = 1 / A(i,i);
    A = S(:,:,i) * A;
         end

        ProdE = 1;
        [~,~,c] = size(E);
        for ii = c:-1:1,
    ProdE = ProdE * E(:,:,ii);
        end
        ProdP = 1;
        [~,~,c] = size(P);
        for ii = c:-1:1,
    ProdP = ProdP * P(:,:,ii);
        end
        ProdS = 1;
        [~,~,d] = size(S);
        for jj = d:-1:1,
    ProdS = ProdS * S(:,:,jj);
        end
        Inverse = ProdS*ProdP*ProdE;
        disp(Inverse);
end


