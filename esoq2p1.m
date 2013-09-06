function [q,loss]=esoq2p1(obs,ref,wt)

% ESOQ2.1 - Attitude EStimate OPtimization (F. L. Markley, 7/15/99).
% Variation of Daniele Mortari's ESOQ2 (Paper AAS 97-167, AAS/AIAA
% Space Flight Mechanics Meeting, Huntsville, AL, February 10-12, 1997),
% with new singularity-avoidance and lambda_max computation logic.
%
% input: obs(3,n) - array of n observation unit vectors
% ref(3,n) - array of n reference unit vectors
% wt(n) - row array of n measurement weights
%
% The columns of obs and ref are assumed to be normalized.
% No assumption is made about the normalization of the weights.
%
% output: q(4) - optimal quaternion
%  loss - optimized value of Wahba's loss function

lam=norm(wt,1);	 % zeroth order approximation to lambda_max

B=[obs(1,:).*wt; obs(2,:).*wt; obs(3,:).*wt]*ref'; trB=trace(B);

% Optimal 180 deg rotation to avoid zero rotation angle singularity
[Bmin,irot]=min([B(1,1) B(2,2) B(3,3) trB]);
switch irot
 case 1, B=[ B(:,1) -B(:,2:3)];			trB=2*Bmin-trB;
 case 2, B=[-B(:,1) B(:,2) -B(:,3)];	trB=2*Bmin-trB;
 case 3, B=[-B(:,1:2) B(:,3)];			trB=2*Bmin-trB;
end

% Compute needed matrices and vectors
S11=2*B(1,1); S23=B(2,3)+B(3,2);
S22=2*B(2,2); S31=B(3,1)+B(1,3);
S33=2*B(3,3); S12=B(1,2)+B(2,1);

z=[B(2,3)-B(3,2); B(3,1)-B(1,3); B(1,2)-B(2,1)];
z12=z(1)*z(1); z22=z(2)*z(2); z32=z(3)*z(3);

% max eigenvalue computation for two observation case
if length(wt)== 2,
   lam0=lam; trB2=trB*trB;
   Sz=[S11 S12 S31; S12 S22 S23; S31 S23 S33]*z;
   aa=trB2-S22*S33+S23*S23-S11*S33+S31*S31-S22*S11+S12*S12;
   bb=trB2+z12+z22+z32; c2=-aa-bb; u=2*sqrt(aa*bb-Sz'*Sz);
   lam=(sqrt(u-c2)+sqrt(-u-c2))/2; loss=lam0-lam;
end
tml=trB-lam; tpl=trB+lam;
M11=tml*(S11-tpl)-z12; M23=tml*S23-z(2)*z(3);
M22=tml*(S22-tpl)-z22; M31=tml*S31-z(3)*z(1);
M33=tml*(S33-tpl)-z32; M12=tml*S12-z(1)*z(2);

if length(wt)~=2,
   m1=[M11 M12 M31]; m2=[M12 M22 M23]; m3=[M31 M23 M33];
   n1=[(S11-2*lam) S12 S31]; n2=[S12 (S22-2*lam) S23]; n3=[S31 S23 (S33-2*lam)];
end

% Compute loss function and rotation axis
e=[M22*M33-M23*M23; M11*M33-M31*M31; M11*M22-M12*M12];
[dummy,k]=max(abs(e));
switch k,
case 1,
   e=[e(1); M31*M23-M12*M33; M12*M23-M31*M22];
   if length(wt)~=2,
      v=cross(m2,n3)'-cross(m3,n2)'; loss=-(m1*e)/(n1*e+m1*v);
      tml=tml+loss; e=e+loss*v;
   end
case 2,
   e=[M31*M23-M12*M33; e(2); M12*M31-M11*M23];
   if length(wt)~=2,
      v=cross(m3,n1)'-cross(m1,n3)'; loss=-(m2*e)/(n2*e+m2*v);
      tml=tml+loss; e=e+loss*v;
   end
case 3,
   e=[M12*M23-M31*M22; M12*M31-M11*M23; e(3)];
   if length(wt)~=2,
      v=cross(m1,n2)'-cross(m2,n1)'; loss=-(m3*e)/(n3*e+m3*v);
      tml=tml+loss; e=e+loss*v;
   end
end

% Quaternion computation in rotated frame
q=[tml*e; -z'*e]; q=q/sqrt(q'*q);

% Undo rotation to get quaternion in input frame
switch irot
 case 1, q=[ q(4);-q(3); q(2);-q(1)];
 case 2, q=[ q(3); q(4);-q(1);-q(2)];
 case 3, q=[-q(2); q(1); q(4);-q(3)];
end
 
