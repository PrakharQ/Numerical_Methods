filename=input('Enter the name of the file from where to take inputs\n','s');
fileid=fopen(filename,'r'); 
formatSpec = '%f';
n = fscanf(fileid,formatSpec,1);
B=fscanf(fileid,formatSpec,[n n]);
max_rel_error=0.001;
A=transpose(B);
x=input('Only the largest eigenvalue (L) or all eigenvalues (A)?\n','s');
if strcmp(x,'L')
    fclose(fileid);
        B=A;
    x=ones(n,1);
    i=0;
    rel_error=100;
    while(abs(rel_error)>max_rel_error)
    x=(B*x)/norm(B*x);
    eigen=(transpose(x)*B*x)/(transpose(x)*x);
    e=eigen;
    if(i>=1)
        rel_error=(e - e_prev)*100/e;
    end
    e_prev=e;
    i=i+1;
    end
    disp(eigen);
    fileID = fopen('outputQ2_Power.txt','w');
    fprintf(fileID,'Eigenvalue\n%f\n\n',eigen);
    fclose(fileID);
    type outputQ2_power.txt
end
if strcmp(x,'A')
    C=A;
    p=1;
    max_no_iteration=50;
    max_rel_error=0.001;
    while( p < max_no_iteration )
       C_prev=C;
    for i = 1:n
       C(:,i) = C(:,i) / norm(C(:,i));
    for j = i+1:1:n
       C(:,j) = C(:,j) - dot(C(:,i),C(:,j))/dot(C(:,i),C(:,i))*C(:,i);
    end
    end

    Q=C;
    R=inv(C)*C_prev;
    C=R*Q;
    p=p+1;
    end
    e=diag(C);
    disp(e);
    fileID = fopen('outputQ2_QR.txt','w');
    fprintf(fileID,'Eigenvalues= \n%f \n \n');
    fprintf(fileID,' %f \n ',e);
    fclose(fileID);
    type outputQ2_QR.txt
end


