a=5;
b=4;
c=5;
d=2;
X0=[0.5;-0.1;0.3;-0.8];
pole1 = -1.5+1.5j;
pole2 = -1.5-1.5j;
pole3 = -2+2j;
pole4 = -2-2j;
X0=[0.5;-0.1;0.3;-0.8];
syms s
phi_s2=(s-pole1)*(s-pole2);
phi_s1=(s-pole3)*(s-pole4);
Ad_lastrow1 = fliplr(-sym2poly(phi_s1));
Ad_lastrow2 = fliplr(-sym2poly(phi_s2));
Ad=[0,1,0,0;
    Ad_lastrow1(1:2),0,0;
    0,0,0,1;
    0,0,Ad_lastrow2(1:2)];
A=[-8.8487+(a-b)/5, -0.0399, -5.55+(c+d)/10, 3.5846;
   -4.5740, 2.5010*((d+5)/(c+5)), -4.3662, -1.1183-(a-c)/20;
   3.7698, 16.1212-c/5, -18.2103+(a+d)/(b+4), 4.4936;
   -8.5645-(a-b)/(c+d+2), 8.3742, -4.4331, -7.7181*(c+5)/(b+5)];
B=[0.0564+b/(10+c), 0.0319;
   0.0165-(c+d-5)/(1000+20*a), -0.02;
   4.4939, 1.5985*(a+10)/(b+12);
   -1.4269, -0.2730];
C=[-3.2988,-2.1932+(10*c+d)/(100+5*a), 0.0370, -0.0109;
   0.2922- a*b/500, -2.150, -0.0104, 0.0163];
R=1;
Q=[50,0,0,0;
   0,50,0,0;
   0,0,50,0;
   0,0,0,50];
Gamma = [A,-B*(inv(R))*(transpose(B)); -Q,-(transpose(A))];
[V,D]=eig(Gamma);
for i= 1:8
    if real(D(i,i)) > 0
        V(:,i)=[0];
    end
end
V(:,all(V==0,1))=[];
v = V(1:4,:);
mu = V(5:8,:);
P=mu*(inv(v));
K=round(((inv(R))*(transpose(B))*P),4);
Atuta = transpose(A);
Btuta = transpose(C);
Wc = [Btuta, Atuta*Btuta,Atuta*Atuta*Btuta,Atuta*Atuta*Atuta*Btuta];
Con=[Wc(:,1),Wc(:,3),Wc(:,2),Wc(:,4)];
    Con_inverse=inv(Con);
    q2T=Con_inverse(2,:);
    q4T=Con_inverse(4,:);
    T=[q2T;
       q2T*Atuta;
       q4T;
       q4T*Atuta];
   Atutabar=T*Atuta*(inv(T));
   Ktutabar=Atutabar-Ad;
   Ktutabar(1,:)=[];
   Ktutabar(2,:)=[];
   Ktuta = Ktutabar*T;
   L = round(transpose(Ktuta),4);

