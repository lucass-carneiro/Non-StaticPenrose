(* ::Package:: *)

(* ::Section:: *)
(*Options (User)*)


(* The coordinate basis *)
ClearAll[coords];
coords = {t, x, y, z};

(* The function to use for simplifications *)
ClearAll[simplifyier];
simplifyier = FullSimplify;


(* ::Section:: *)
(*Spacetime metric (User)*)


funcs = {t[bhIndex][T,X,Y,Z], x[bhIndex][T,X,Y,Z], y[bhIndex][T,X,Y,Z], z[bhIndex][T,X,Y,Z]}
J = FullSimplify[Table[D[funcs[[i]],{T,X,Y,Z}[[j]]], {i, 1, 3}, {j, 1, 3}]];
iJ = FullSimplify[Inverse[J]]


(* Auxiliary metric functions *)
ClearAll[l];
l[bhIndex_][x_,y_,z_]:= {1, l1[bhIndex][x,y,z], l2[bhIndex][x,y,z], l3[bhIndex][x,y,z] };
\[ScriptCapitalH][bhIndex_][x_,y_,z_]:= Table[2 * H[bhIndex] * l[bhIndex][x,y,z][[a]]l[bhIndex][x,y,z][[b]], {a, 1, 3}, {b, 1,3}]

\[ScriptCapitalH][1][x,y,z]//MatrixForm


eta = DiagonalMatrix[{-1,1,1,1}];

(* A 4x4 matrix with the metric to generate *)
ClearAll[ll4metric];
ll4metric =FullSimplify[Table[eta[[mu,nu]] + 2*H[x,y,z] * l[[mu]] * l[[nu]], {mu,1,4}, {nu,1,4}]];

$Assumptions = M > 0 &&
               a > 0 &&
               r[x,y,z] >= 0;


(* ::Section::Closed:: *)
(*ADM components*)


(* Assume that the coordinates are all real *)
$Assumptions = $Assumptions && And@@(Element[#,Reals]&/@coords);

(* The inverse metric. The user might want to set this directly for speed *)
ClearAll[uu4metric];
uu4metric = simplifyier[Inverse[ll4metric]];

(* Spatial coordinates *)
ClearAll[spaceCoords];
spaceCoords = coords[[2;;4]];

ClearAll[llsmetric,uusmetric];
llsmetric = simplifyier[ll4metric[[2;;4, 2;;4]]];
uusmetric = simplifyier[Inverse[llsmetric]];

ClearAll[lapse];
lapse = simplifyier[Sqrt[-1/uu4metric[[1,1]]]];

ClearAll[lshift,ushift];
lshift = simplifyier[ll4metric[[1, 2;;4]]];
ushift = simplifyier[lapse^2 * uu4metric[[1, 2;;4]]];

ClearAll[\[CapitalGamma],Dd];
\[CapitalGamma][i_,j_,k_] := simplifyier[Sum[1/2 * uusmetric[[i,l]] * ( D[llsmetric[[l,j]],spaceCoords[[k]]] + D[llsmetric[[l,k]],spaceCoords[[j]]] - D[llsmetric[[j,k]],spaceCoords[[l]]] ),{l,1,3}]];
Dd[b_,c_,f_] := simplifyier[D[f[[c]], spaceCoords[[b]]] - Sum[\[CapitalGamma][d,b,c]f[[d]],{d,1,3}]];

ClearAll[llextrinsic,ulextrinsic];
llextrinsic = simplifyier[1/(2*lapse)*Table[-D[llsmetric[[i,j]],coords[[1]]] + Dd[i,j,lshift] + Dd[j,i,lshift], {i,1,3}, {j, 1,3}]];
ulextrinsic = simplifyier[Table[Sum[uusmetric[[i,k]]*llextrinsic[[k,j]],{k,1,3}],{i,1,3},{j,1,3}]];


(* ::Section::Closed:: *)
(*ADM Derivatives*)


ClearAll[gradlapse];
gradlapse = simplifyier[Table[D[lapse,spaceCoords[[i]]], {i,1,3}]];

ClearAll[gradushift];
gradushift = Table[D[ushift[[j]],spaceCoords[[i]]],{i,1,3},{j,1,3}];


(* ::Section::Closed:: *)
(*Printout*)


ClearAll[rules];
rules={
  "H(x,y,z)"->"H",
  
  "l1(x,y,z)"->"l1",
  "l2(x,y,z)"->"l2",
  "l3(x,y,z)"->"l3",
  
  "r(x,y,z)" -> "r",
  
  "Derivative(1,0,0)(H)(x,y,z)" -> "dH_dx",
  "Derivative(0,1,0)(H)(x,y,z)" -> "dH_dy",
  "Derivative(0,0,1)(H)(x,y,z)" -> "dH_dz",
  
  "Derivative(1,0,0)(l1)(x,y,z)" -> "dl1_dx",
  "Derivative(0,1,0)(l1)(x,y,z)" -> "dl1_dy",
  "Derivative(0,0,1)(l1)(x,y,z)" -> "dl1_dz",
  
  "Derivative(1,0,0)(l2)(x,y,z)" -> "dl2_dx",
  "Derivative(0,1,0)(l2)(x,y,z)" -> "dl2_dy",
  "Derivative(0,0,1)(l2)(x,y,z)" -> "dl2_dz",
  
  "Derivative(1,0,0)(l3)(x,y,z)" -> "dl3_dx",
  "Derivative(0,1,0)(l3)(x,y,z)" -> "dl3_dy",
  "Derivative(0,0,1)(l3)(x,y,z)" -> "dl3_dz",
  
  "Derivative(1,0,0)(r)(x,y,z)" -> "dr_dx",
  "Derivative(0,1,0)(r)(x,y,z)" -> "dr_dy",
  "Derivative(0,0,1)(r)(x,y,z)" -> "dr_dz"
};


StringReplace[ToString[FullSimplify[lapse],CForm],rules]<>";"


StringReplace[ToString[FullSimplify[lshift[[1]]],CForm],rules]<>";"
StringReplace[ToString[FullSimplify[lshift[[2]]],CForm],rules]<>";"
StringReplace[ToString[FullSimplify[lshift[[3]]],CForm],rules]<>";"


StringReplace[ToString[FullSimplify[ushift[[1]]],CForm],rules]<>";"
StringReplace[ToString[FullSimplify[ushift[[2]]],CForm],rules]<>";"
StringReplace[ToString[FullSimplify[ushift[[3]]],CForm],rules]<>";"


StringReplace[ToString[FullSimplify[llextrinsic[[1,1]]],CForm],rules]<>";"
StringReplace[ToString[FullSimplify[llextrinsic[[1,2]]],CForm],rules]<>";"
StringReplace[ToString[FullSimplify[llextrinsic[[1,3]]],CForm],rules]<>";"
StringReplace[ToString[FullSimplify[llextrinsic[[2,2]]],CForm],rules]<>";"
StringReplace[ToString[FullSimplify[llextrinsic[[2,3]]],CForm],rules]<>";"
StringReplace[ToString[FullSimplify[llextrinsic[[3,3]]],CForm],rules]<>";"


StringReplace[ToString[Simplify[ulextrinsic[[1,1]]],CForm],rules]<>";"
StringReplace[ToString[Simplify[ulextrinsic[[1,2]]],CForm],rules]<>";"
StringReplace[ToString[Simplify[ulextrinsic[[1,3]]],CForm],rules]<>";"
StringReplace[ToString[Simplify[ulextrinsic[[2,1]]],CForm],rules]<>";"
StringReplace[ToString[Simplify[ulextrinsic[[2,2]]],CForm],rules]<>";"
StringReplace[ToString[Simplify[ulextrinsic[[2,3]]],CForm],rules]<>";"
StringReplace[ToString[Simplify[ulextrinsic[[3,1]]],CForm],rules]<>";"
StringReplace[ToString[Simplify[ulextrinsic[[3,2]]],CForm],rules]<>";"
StringReplace[ToString[Simplify[ulextrinsic[[3,3]]],CForm],rules]<>";"


StringReplace[ToString[FullSimplify[gradlapse[[1]]],CForm],rules]<>";"
StringReplace[ToString[FullSimplify[gradlapse[[2]]],CForm],rules]<>";"
StringReplace[ToString[FullSimplify[gradlapse[[3]]],CForm],rules]<>";"


StringReplace[ToString[FullSimplify[gradushift[[1,1]]],CForm],rules]<>";"
StringReplace[ToString[FullSimplify[gradushift[[1,2]]],CForm],rules]<>";"
StringReplace[ToString[FullSimplify[gradushift[[1,3]]],CForm],rules]<>";"
StringReplace[ToString[FullSimplify[gradushift[[2,1]]],CForm],rules]<>";"
StringReplace[ToString[FullSimplify[gradushift[[2,2]]],CForm],rules]<>";"
StringReplace[ToString[FullSimplify[gradushift[[2,3]]],CForm],rules]<>";"
StringReplace[ToString[FullSimplify[gradushift[[3,1]]],CForm],rules]<>";"
StringReplace[ToString[FullSimplify[gradushift[[3,2]]],CForm],rules]<>";"
StringReplace[ToString[FullSimplify[gradushift[[3,3]]],CForm],rules]<>";"


StringReplace[ToString[FullSimplify[llsmetric[[1,1]]],CForm],rules]<>";"
StringReplace[ToString[FullSimplify[llsmetric[[1,2]]],CForm],rules]<>";"
StringReplace[ToString[FullSimplify[llsmetric[[1,3]]],CForm],rules]<>";"
StringReplace[ToString[FullSimplify[llsmetric[[2,2]]],CForm],rules]<>";"
StringReplace[ToString[FullSimplify[llsmetric[[2,3]]],CForm],rules]<>";"
StringReplace[ToString[FullSimplify[llsmetric[[3,3]]],CForm],rules]<>";"


StringReplace[ToString[FullSimplify[uusmetric[[1,1]]],CForm],rules]<>";"
StringReplace[ToString[FullSimplify[uusmetric[[1,2]]],CForm],rules]<>";"
StringReplace[ToString[FullSimplify[uusmetric[[1,3]]],CForm],rules]<>";"
StringReplace[ToString[FullSimplify[uusmetric[[2,2]]],CForm],rules]<>";"
StringReplace[ToString[FullSimplify[uusmetric[[2,3]]],CForm],rules]<>";"
StringReplace[ToString[FullSimplify[uusmetric[[3,3]]],CForm],rules]<>";"


StringReplace[ToString[FullSimplify[\[CapitalGamma][1,1,1]],CForm],rules]<>";"
StringReplace[ToString[FullSimplify[\[CapitalGamma][1,1,2]],CForm],rules]<>";"
StringReplace[ToString[FullSimplify[\[CapitalGamma][1,1,3]],CForm],rules]<>";"
StringReplace[ToString[FullSimplify[\[CapitalGamma][1,2,2]],CForm],rules]<>";"
StringReplace[ToString[FullSimplify[\[CapitalGamma][1,2,3]],CForm],rules]<>";"
StringReplace[ToString[FullSimplify[\[CapitalGamma][1,3,3]],CForm],rules]<>";"

StringReplace[ToString[FullSimplify[\[CapitalGamma][2,1,1]],CForm],rules]<>";"
StringReplace[ToString[FullSimplify[\[CapitalGamma][2,1,2]],CForm],rules]<>";"
StringReplace[ToString[FullSimplify[\[CapitalGamma][2,1,3]],CForm],rules]<>";"
StringReplace[ToString[FullSimplify[\[CapitalGamma][2,2,2]],CForm],rules]<>";"
StringReplace[ToString[FullSimplify[\[CapitalGamma][2,2,3]],CForm],rules]<>";"
StringReplace[ToString[FullSimplify[\[CapitalGamma][2,3,3]],CForm],rules]<>";"

StringReplace[ToString[FullSimplify[\[CapitalGamma][3,1,1]],CForm],rules]<>";"
StringReplace[ToString[FullSimplify[\[CapitalGamma][3,1,2]],CForm],rules]<>";"
StringReplace[ToString[FullSimplify[\[CapitalGamma][3,1,3]],CForm],rules]<>";"
StringReplace[ToString[FullSimplify[\[CapitalGamma][3,2,2]],CForm],rules]<>";"
StringReplace[ToString[FullSimplify[\[CapitalGamma][3,2,3]],CForm],rules]<>";"
StringReplace[ToString[FullSimplify[\[CapitalGamma][3,3,3]],CForm],rules]<>";"


(* ::Section::Closed:: *)
(*Auxiliary functions Printouts*)


ClearAll[r]
r[x_,y_,z_]=rvar/.Solve[(x^2+y^2)/(rvar^2+a^2)+z^2/rvar^2==1,rvar][[4]]//FullSimplify;

ToString[FullSimplify[r[x,y,z]],CForm]
ToString[FullSimplify[D[r[x,y,z],x]],CForm]
ToString[FullSimplify[D[r[x,y,z],y]],CForm]
ToString[FullSimplify[D[r[x,y,z],z]],CForm]

ClearAll[r]


ClearAll[H];
H[x_,y_,z_]:=(M*r[x,y,z])/(r[x,y,z]^2+a^2*(z/r[x,y,z])^2)

StringReplace[ToString[FullSimplify[H[x,y,z]],CForm],rules]
StringReplace[ToString[FullSimplify[D[H[x,y,z],x]],CForm],rules]
StringReplace[ToString[FullSimplify[D[H[x,y,z],y]],CForm],rules]
StringReplace[ToString[FullSimplify[D[H[x,y,z],z]],CForm],rules]

ClearAll[H];


ClearAll[l1];
l1[x_,y_,z_]:=(r[x,y,z]*x+a*y)/(r[x,y,z]^2+a^2)

StringReplace[ToString[FullSimplify[l1[x,y,z]],CForm],rules]
StringReplace[ToString[FullSimplify[D[l1[x,y,z],x]],CForm],rules]
StringReplace[ToString[FullSimplify[D[l1[x,y,z],y]],CForm],rules]
StringReplace[ToString[FullSimplify[D[l1[x,y,z],z]],CForm],rules]

ClearAll[l1];


ClearAll[l2];
l2[x_,y_,z_]:=(r[x,y,z]*y-a*x)/(r[x,y,z]^2+a^2)

StringReplace[ToString[FullSimplify[l2[x,y,z]],CForm],rules]
StringReplace[ToString[FullSimplify[D[l2[x,y,z],x]],CForm],rules]
StringReplace[ToString[FullSimplify[D[l2[x,y,z],y]],CForm],rules]
StringReplace[ToString[FullSimplify[D[l2[x,y,z],z]],CForm],rules]

ClearAll[l2];


ClearAll[l3];
l3[x_,y_,z_]:=z/r[x,y,z]

StringReplace[ToString[FullSimplify[l3[x,y,z]],CForm],rules]
StringReplace[ToString[FullSimplify[D[l3[x,y,z],x]],CForm],rules]
StringReplace[ToString[FullSimplify[D[l3[x,y,z],y]],CForm],rules]
StringReplace[ToString[FullSimplify[D[l3[x,y,z],z]],CForm],rules]

ClearAll[l3];





gradlapse[[2]]//FullSimplify


(* ::Section:: *)
(*Pointwise numeric printouts*)


ClearAll[r];
ClearAll[H];
ClearAll[l1];
ClearAll[l2];
ClearAll[l3];

up=En*{1/lapse,V0-ushift[[1]]/lapse,V1-ushift[[2]]/lapse,V2-ushift[[3]]/lapse};
\[Xi]={1,0,0,0};
Eg=FullSimplify[-Sum[ll4metric[[a,b]]up[[a]]\[Xi][[b]],{a,1,4},{b,1,4}]];

(*r[x_,y_,z_]=rvar/.Solve[(x^2+y^2)/(rvar^2+a^2)+z^2/rvar^2==1,rvar][[4]]//FullSimplify;
H[x_,y_,z_]:=(M*r[x,y,z])/(r[x,y,z]^2+a^2*(z/r[x,y,z])^2)
l1[x_,y_,z_]:=(r[x,y,z]*x+a*y)/(r[x,y,z]^2+a^2)
l2[x_,y_,z_]:=(r[x,y,z]*y-a*x)/(r[x,y,z]^2+a^2)
l3[x_,y_,z_]:=z/r[x,y,z]*)


ClearAll[r];
ClearAll[H];
ClearAll[l1];
ClearAll[l2];
ClearAll[l3];


Eg=En ((2 M (a^2 (V0 x+V1 y-V2 z)+Sqrt[2] a (V1 x-V0 y) Sqrt[-a^2+x^2+y^2+z^2+Sqrt[(-a^2+x^2+y^2)^2+2 (a^2+x^2+y^2) z^2+z^4]]-(V0 x+V1 y+V2 z) (x^2+y^2+z^2+Sqrt[4 a^2 z^2+(-a^2+x^2+y^2+z^2)^2])))/(Sqrt[4 a^2 z^2+(-a^2+x^2+y^2+z^2)^2] (a^2+x^2+y^2+z^2+Sqrt[4 a^2 z^2+(-a^2+x^2+y^2+z^2)^2]))+1/Sqrt[1+Sqrt[2] M Sqrt[(-a^2+x^2+y^2+z^2+Sqrt[4 a^2 z^2+(-a^2+x^2+y^2+z^2)^2])/(4 a^2 z^2+(-a^2+x^2+y^2+z^2)^2)]]);


Reduce[FullSimplify[Eg//.{M->1,a->9/10,y->0,z->0,x->-19/10,En->1/2,V0->-1/2}]==-10^-2,Reals]//N
