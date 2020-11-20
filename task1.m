a=5;
b=4;
c=5;
d=2;
pole1 = -0.25+0.25j;
pole2 = -0.25-0.25j;
pole3 = -1+1j;
pole4 = -1-1j;
X0=[0.5;-0.1;0.3;-0.8];
syms s
phi_s2=(s-pole1)*(s-pole2);
phi_s1=(s-pole3)*(s-pole4);
Ad_lastrow1 = fliplr(-sym2poly(phi_s1));
Ad_lastrow2 = fliplr(-sym2poly(phi_s2));
phi_s = (s-pole1)*(s-pole2)*(s-pole3)*(s-pole4);
Ad_lastrow = fliplr(-sym2poly(phi_s));
Ad=[0,1,0,0;
    Ad_lastrow1(1:2),0,0;
    0,0,0,1;
    0,0,Ad_lastrow2(1:2)];
%Ad=[0,1,0,0;
 %   0,0,1,0;
  %  0,0,0,1;
   % Ad_lastrow(1:4)];
A=[-8.8487+(a-b)/5, -0.0399, -5.55+(c+d)/10, 3.5846;
   -4.5740, 2.5010*((d+5)/(c+5)), -4.3662, -1.1183-(a-c)/20;
   3.7698, 16.1212-c/5, -18.2103+(a+d)/(b+4), 4.4936;
   -8.5645-(a-b)/(c+d+2), 8.3742, -4.4331, -7.7181*(c+5)/(b+5)];
B=[0.0564+b/(10+c), 0.0319;
   0.0165-(c+d-5)/(1000+20*a), -0.02;
   4.4939, 1.5985*(a+10)/(b+12);
   -1.4269, -0.2730];
C=[-3.2988,-2.1932+(10*c+d)/(100+5*a), 0.0370, -0.0109;
   0.2922- a*b/500, -7.1506, -0.0104, 0.0163];

Wc=[B,A*B,A*A*B,A*A*A*B];
if rank(Wc)~= 4 
    print("not controllable")
else
    Con=[Wc(:,1),Wc(:,3),Wc(:,2),Wc(:,4)];
    Con_inverse=inv(Con);
    q2T=Con_inverse(2,:);
    q4T=Con_inverse(4,:);
    T=[q2T;
       q2T*A;
       q4T;
       q4T*A];
   Abar=T*A*(inv(T));
   Kbar=Abar-Ad;
   Kbar(1,:)=[];
   Kbar(2,:)=[];
   K=Kbar*T;
   Bbar=T*B;
end
    
