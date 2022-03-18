(* ::Package:: *)

(* ::Section:: *)
(*Options (User)*)


(* The coordinate basis *)
ClearAll[coords];
coords = {t, x, y, z};

(* The function to use for simplifications *)
ClearAll[simplifyier];
simplifyier = Simplify;


(* ::Section:: *)
(*Spacetime metric (User)*)


(* Auxiliary metric functions *)
ClearAll[f,k,eta];
f[x,y,z] := (2 * M * r[x,y,z]^3)/(r[x,y,z]^4 + a^2 * z^2);
k = {1, (r[x,y,z] * x + a * y)/(r[x,y,z]^2 + a^2), (r[x,y,z] * y - a * x)/(r[x,y,z]^2 + a^2), z/r[x,y,z]};
eta = DiagonalMatrix[{-1,1,1,1}];

(* A 4x4 matrix with the metric to generate *)
ClearAll[ll4metric];
ll4metric = Table[eta[[mu,nu]] + f[x,y,z] * k[[mu]] * k[[nu]], {mu,1,4}, {nu,1,4}];

$Assumptions = M > 0 && a > 0 && r[x,y,z] >= 0;


ClearAll[r];
(*r[x_,y_,z_]:=Sqrt[Sqrt[(-a^2+x^2+y^2+z^2)^2+4*a^2*z^2]-a^2+x^2+y^2+z^2]/Sqrt[2];*)
ClearAll[r];


(* ::Section:: *)
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
\[CapitalGamma][a_,b_,c_] := simplifyier[1/2 * Sum[uusmetric[[a,d]] * ( D[llsmetric[[d,b]],spaceCoords[[c]]] + D[llsmetric[[d,c]],spaceCoords[[b]]] - D[llsmetric[[b,c]],spaceCoords[[d]]] ),{d,1,3}]];
Dd[b_,c_,f_] := simplifyier[D[f[[c]], spaceCoords[[b]]] - Sum[\[CapitalGamma][d,b,c]f[[d]],{d,1,3}]];

ClearAll[llextrinsic,ulextrinsic];
llextrinsic = simplifyier[1/(2*lapse)*Table[-D[llsmetric[[i,j]],coords[[1]]] + Dd[i,j,lshift] + Dd[j,i,lshift], {i,1,3}, {j, 1,3}]];
ulextrinsic = simplifyier[Table[Sum[uusmetric[[i,k]]*llextrinsic[[k,j]],{k,1,3}],{i,1,3},{j,1,3}]];


(* ::Section:: *)
(*ADM Derivatives*)


ClearAll[gradlapse];
gradlapse = simplifyier[Table[D[lapse,spaceCoords[[i]]], {i,1,3}]];

ClearAll[gradushift];
gradushift = Table[D[ushift[[i]],spaceCoords[[j]]],{i,1,3},{j,1,3}];


(* ::Section:: *)
(*Printout*)


ClearAll[rules];
rules={
  "r(x,y,z)" -> "r",
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


StringReplace[ToString[FullSimplify[(Sqrt[Sqrt[(-a^2+x^2+y^2+z^2)^2+4*a^2*z^2]-a^2+x^2+y^2+z^2]/Sqrt[2])],CForm],rules]<>";"


StringReplace[ToString[FullSimplify[D[(Sqrt[Sqrt[(-a^2+x^2+y^2+z^2)^2+4*a^2*z^2]-a^2+x^2+y^2+z^2]/Sqrt[2]),x]],CForm],rules]<>";"


StringReplace[ToString[FullSimplify[D[(Sqrt[Sqrt[(-a^2+x^2+y^2+z^2)^2+4*a^2*z^2]-a^2+x^2+y^2+z^2]/Sqrt[2]),y]],CForm],rules]<>";"


StringReplace[ToString[FullSimplify[D[(Sqrt[Sqrt[(-a^2+x^2+y^2+z^2)^2+4*a^2*z^2]-a^2+x^2+y^2+z^2]/Sqrt[2]),z]],CForm],rules]<>";"
