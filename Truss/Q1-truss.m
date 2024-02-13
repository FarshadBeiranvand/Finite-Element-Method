clc;
clear all;


%readings datas from truss input.xlsx

elnod=xlsread('truss input.xlsx','element nodes');
f=xlsread('truss input.xlsx','force');
ncor=xlsread('truss input.xlsx','nodecoordinates');
dt=xlsread('truss input.xlsx','deltat');
A=xlsread('truss input.xlsx','Area');
E=xlsread('truss input.xlsx','Elasticity');
alfa=xlsread('truss input.xlsx','alfa');
fd=xlsread('truss input.xlsx','fix displacements');
intd=xlsread('truss input.xlsx','initial displacements');
is=xlsread('truss input.xlsx','inclined support');

%nel= number of elements, nnod= number of nodes
v=size(elnod);
nel=v(1,1);
nnod=max(max(elnod(:,2:3)));

%calculate teta=t and length=l of each element
l=zeros(nel,1);
t=zeros(nel,1);
for n=1:nel
    i=elnod(n,2);
    j=elnod(n,3);
    l(n,1)=((ncor(i,1)-ncor(j,1))^2+(ncor(i,2)-ncor(j,2))^2)^0.5;
    t(n,1)=atan ((ncor(j,2)-ncor(i,2))/(ncor(j,1)-ncor(i,1)));
end

%conectivity table
for i=1:nel
    con(i,:)=[2*elnod(i,2)-1 2*elnod(i,2) 2*elnod(i,3)-1 2*elnod(i,3)];
end

%matrix of each element
for k=1:nel
    s=sin(t(k,1));
    c=cos(t(k,1));
    ael=A(k)*E(k)/l(k);
    K(:,:,k)=ael.*[c^2 c*s -c^2 -c*s;c*s s^2 -c*s -s^2;-c^2 -c*s c^2 c*s;c*s -s^2 c*s s^2];
end

%thermal force of each element
for i=1:nel
    fff=A(1,i)*E(1,i)*alfa*dt(1,i);
    fth(:,:,i)=fff.*[-cos(t(i));-sin(t(i));cos(t(i));sin(t(i))];
end

%assembly of thermal force

fthg=zeros(nnod*2,1);
for k=1:nel
    for i=1:4
        I=con(k,i);
        fthg(I,1)=fthg(I,1)+fth(i,1,k);
    end
end

%assembly of k matrix
kg=zeros(nnod*2,nnod*2);
for k=1:nel
    for i=1:4
        for j=1:4
            I=con(k,i);
            J=con(k,j);
            kg(I,J)=kg(I,J)+K(i,j,k);
        end
    end
end

%inclined support
[a1 b1]=size(is);
T=eye(2*nnod,2*nnod);
for i=1:a1
    z=is(i,2)*pi/180;
    T(2*is(i,1)-1,2*is(i,1)-1)=cos(z);
    T(2*is(i,1)-1,2*is(i,1))=sin(z);
    T(2*is(i,1),2*is(i,1)-1)=-sin(z);
    T(2*is(i,1),2*is(i,1))=cos(z);
end

kg=T*kg*(T');
fthg=(T')*fthg;
KG=kg;

FG=fthg+f;
FF=FG;
for i=1:length(fd)
    kg(fd(i),:)=0;
    kg(:,fd(i))=0;
    kg(fd(i),fd(i))=1;
    FG(fd(i),:)=0;
end
 
%initial displacement
[a2 b2]=size(intd);
for i=1:a2
    FG=FG-intd(i,2)*kg(:,intd(i,1));
    FG(intd(i,1),1)=intd(i,2);
    kg(intd(i,1),:)=0;
    kg(:,intd(i,1))=0;
    kg(intd(i,1),intd(i,1))=1;
end

%calculate displacements

d=inv(kg)*FG;

F=KG*d;

%supports fprce=R
R=KG*d-FF;

%stress of each element=st
del=zeros(4,1);
TT=zeros(1,4);
sigma=zeros(nel,1);
for n=1:nel
    i=elnod(n,2);
    j=elnod(n,3);
    c=cos(t(n,1));
    s=sin(t(n,1));
    TT=[-c -s c s];
    del=[d(2*i-1,1);d(2*i,1);d(2*j-1,1);d(2*j,1)];
    sigma(n,1)=(E(1,n)/l(n,1))*TT*del;
end
%force for each element
fel=zeros(nel,1);
for i=1:nel
    fel(i,1)=sigma(i,1)*A(1,i);
end


%increasing of element length
deltal=zeros(nel,1);
for n=1:nel
    i=elnod(n,2);
    j=elnod(n,3);
    deltal(n,1)=sqrt((d(2*i-1)-d(2*j-1))^2+(d(2*i)-d(2*j))^2);
end


%reshape R
rr=zeros(nnod,3);
for i=1:nnod
    rr(i,1)=i;
    rr(i,2)=R(2*i-1,1);
    rr(i,3)=R(2*i,1);
end

%reshape d
dd=zeros(nnod,3);
for i=1:nnod
    dd(i,1)=i;
    dd(i,2)=d(2*i-1,1);
    dd(i,3)=d(2*i,1);
end





%writing the results in excel


xlswrite ('truss output.xlsx',dd,'displacement','A2');
xlswrite ('truss output.xlsx',rr,'support forces','A2');
xlswrite ('truss output.xlsx',fel,'force for each element','A1');    
xlswrite ('truss output.xlsx',deltal,'delta length','A1');






