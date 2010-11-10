%
% M-file for matlab.
%
% Calculates the regular icosahedron vertices for the TERRA grid.
%

a=1.;
tau=0.5*(sqrt(5)+1);
rho=tau-1;
u=a/(sqrt(1+rho^2));
v=rho*a/(sqrt(1+rho^2));
%A=[
% u  v  0,
% u -v  0,
%-u  v  0,
%-u -v  0,
% 0  u  v,
% 0  u -v,
% 0 -u  v,
% 0 -u -v,
% v  0  u,
% v  0 -u,
%-v  0  u,
%-v  0 -u];

A=[
 v  0  u,
 u  v  0,
 0  u  v,
-v  0  u,
 0 -u  v,
 u -v  0,
 v  0 -u,
 0  u -v,
-u  v  0,
-u -v  0,
 0 -u -v,
-v  0 -u];



Beta=atan(v/u);

Ry=[
    cos(Beta)   0  -sin(Beta),
       0        1      0     ,
    sin(Beta)  0    cos(Beta)];

for i=1:12
	  Ad(i,:)=Ry*transpose(A(i,:));
end

A

Ad

scatter3(A(:,1),A(:,2),A(:,3))
hold on
scatter3(Ad(:,1),Ad(:,2),Ad(:,3),'filled')
xlabel('x')
ylabel('y')
zlabel('z')
