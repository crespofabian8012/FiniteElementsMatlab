load('coordenadas.mat')

load('triangulos.mat')

load('frontera.mat')

%%%elementos y triangulos

figure
plot(C(:,1), C(:,2),'bo')
title('Puntos de discretización')
figure
for j= 1:size(E,1)
    hold on
    plot(C(E(j,:),1),C(E(j,:),2),'Color',[.6 0 0])
    hold on
    
end 
title('Discretización de triángulos')

pointsFrontier=C(unique(F),:)
figure
plot(pointsFrontier(:,1), pointsFrontier(:,2),'bo')
title('Puntos de frontera')


FreeNodes= setdiff(1: size(C,1), unique(F)) 
interiorPoints=C(FreeNodes,:)
figure
plot(interiorPoints(:,1), interiorPoints(:,2),'bo')
title('Puntos interiores')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MatrizMasa=sparse(size(C,1),size(C,1));
A = sparse(size(C,1),size(C,1));
%partes derechas
f1 = sparse(size(C,1),1);
f2 = sparse(size(C,1),1);

%matriz rigidez o stiffness
for j=1:size(E,1)
   A(E(j,:),E(j,:))=A(E(j,:),E(j,:))+stima3(C(E(j,:),:));
end

%matriz masa
for j=1:size(E,1)
   MatrizMasa(E(j,:),E(j,:))=MatrizMasa(E(j,:),E(j,:))+det([1,1,1; C(E(j,:),:)'])*[2,1,1;1,2,1;1,1,2]/24;
end

r=det([1,1,1; C(E(1,:),:)']);
z1= inline('sin(pi*x)*sin(pi*y)*(2*pi^2+1)','x','y');
z2=inline('sin(2*pi*x)*cos(2*pi*y)*(8*pi^2+1)','x','y');
g=inline('sin(2*pi*x)*cos(2*pi*y)','x','y');
for j = 1:size(E,1)
 v=sum(C(E(j,:),:))/3;
 f1(E(j,:)) = f1(E(j,:))+r*z1(v(1),v(2))/6;
 f2(E(j,:)) = f2(E(j,:))+r*z2(v(1),v(2))/6 ;
end
 
A(unique(F),:)=sparse(length(A(unique(F),unique(F))),size(C,1));
A(unique(F),unique(F))=speye(length(A(unique(F),unique(F))));
 
f1(unique(F))=sparse(length(f1(unique(F))),1);

f2(unique(F))=sparse(length(f2(unique(F))),1);

pointsFrontier=C(unique(F),:);

u= sparse(size(C,1),1);
u(unique(F))=u_d(C(unique(F),:));
f2=f2-(A+MatrizMasa)*u;
  
u2(FreeNodes)=(A(FreeNodes, FreeNodes)+MatrizMasa(FreeNodes, FreeNodes))\f2(FreeNodes);
%u2(FreeNodes)= u2(FreeNodes)+u_d(C(FreeNodes,:));
u2(unique(F))=u_d(C(unique(F),:));



u1=(A+MatrizMasa)\f1;
figure
trisurf(E,C(:,1),C(:,2),full(u1'));
title('Solucion FEM  problema 1');

 figure
 trisurf(E,C(:,1),C(:,2),full(u2'));
 title('Solucion FEM  problema 2');
%problema 1
ut1=@(x,y)(sin(pi.*x).*sin(pi.*y));

%problema 2
ut2=@(x,y)(sin(2*pi.*x).*cos(2*pi.*y));
 h=0.01;
 [x1,y1]=meshgrid(h:h:1-h,h:h:1-h);
 UF=ut1(C(:,1),C(:,2));
 Ufe=ut1(x1,y1);
 figure
 surf(Ufe);
title('Solucion analitica problema 1');
display('Para el problema 1');
display('Con L2');
 errorL2=(u1-UF)'*MatrizMasa*(u1-UF);
 display(errorL2);
 display('Número de condición de matriz de  masa');
display(condest(MatrizMasa));
display('Con H1');
errorH1=(u1-UF)'*A*(u1-UF);
display(errorH1);
error=(u1-UF)'*(A+MatrizMasa)*(u1-UF);
 display('Número de condición de matriz de  rigidez');
display(condest(A));
display('Número de condición de suma de matrices');
display(condest(MatrizMasa+A));


UG=ut2(C(:,1),C(:,2));
 Uge=ut2(x1,y1);
figure
 surf(Uge);
title('Solucion analitica problema 2');

display('Para el problema 2');
display('Con L2');
errorL2=(u2'-UG)'*(MatrizMasa)*(u2'-UG);
 display(errorL2);
 display('Número de condicion de matriz de masa');
display(condest(MatrizMasa));
display('Con H1');
errorH1=(u2'-UG)'*A*(u2'-UG);
 display(errorH1);
 error=(u2'-UG)'*(A+MatrizMasa)*(u2'-UG);
 display(error);
 display('Número de condicion de matriz de  rigidez');
display(condest(A));
display('Número de condicion de suma de matrices');
display(condest(MatrizMasa+A));


