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


(* A 4x4 matrix with the metric to generate *)
ClearAll[g];
g = DiagonalMatrix[{-(1 - 2M/(4r[x,y,z]))^2 * (1 + 2*M/(4r[x,y,z]))^(-2), (1 + 2*M/(4r[x,y,z]))^4, (1 + 2*M/(4r[x,y,z]))^4, (1 + 2*M/(4r[x,y,z]))^4}]

$Assumptions = M > 0 && r[x,y,z] > 2 * M/4;


(* ::Section::Closed:: *)
(*ADM components*)


(* Assume that the coordinates are all real *)
$Assumptions = $Assumptions && And@@(Element[#,Reals]&/@coords);

(* The inverse metric. The user might want to set this directly for speed *)
ClearAll[ig];
ig = Inverse[g];

(* Spatial coordinates *)
ClearAll[spaceCoords];
spaceCoords = coords[[2;;4]];

ClearAll[\[Gamma],i\[Gamma]];
\[Gamma] = g[[2;;4, 2;;4]]//simplifyier;
i\[Gamma] = Inverse[\[Gamma]]//simplifyier;

ClearAll[\[Alpha]];
\[Alpha] = Sqrt[-1/ig[[1,1]]]//simplifyier;

ClearAll[\[Beta],i\[Beta]];
\[Beta] = g[[1, 2;;4]]//simplifyier;
i\[Beta] = ig[[1, 2;;4]]/-ig[[1,1]]//simplifyier;

ClearAll[\[CapitalGamma],Dd];
\[CapitalGamma][a_,b_,c_]:=Sum[1/2 * i\[Gamma][[a,d]] * (D[\[Gamma][[d, b]], spaceCoords[[c]]] + D[\[Gamma][[d, c]], spaceCoords[[b]]] - D[\[Gamma][[b, c]], spaceCoords[[d]]]), {d,1,3}]//simplifyier;
Dd[i_,j_,f_]:=D[f[[j]], spaceCoords[[i]]] - Sum[\[CapitalGamma][k,i,j]f[[k]], {k,1,3}]//simplifyier;

ClearAll[K,mixedK];
K = 1/(2*\[Alpha])*Table[-D[\[Gamma][[i,j]],coords[[1]]] + Dd[i,j,\[Beta]] + Dd[j,i,\[Beta]], {i,1,3}, {j, 1,3}]//simplifyier;
mixedK = Table[Sum[i\[Gamma][[i,k]]*K[[k,j]],{k,1,3}],{i,1,3},{j,1,3}];


(* ::Section::Closed:: *)
(*ADM Derivatives*)


ClearAll[grad\[Alpha]];
grad\[Alpha] = Table[D[\[Alpha],spaceCoords[[i]]], {i,1,3}]//simplifyier;

ClearAll[gradi\[Beta]1,gradi\[Beta]2,gradi\[Beta]3];
gradi\[Beta]1 = Table[D[i\[Beta][[1]],spaceCoords[[i]]], {i,1,3}]//simplifyier;
gradi\[Beta]2 = Table[D[i\[Beta][[2]],spaceCoords[[i]]], {i,1,3}]//simplifyier;
gradi\[Beta]3 = Table[D[i\[Beta][[3]],spaceCoords[[i]]], {i,1,3}]//simplifyier;


(* ::Section:: *)
(*Functions*)


ClearAll[r];
r[x_,y_,z_]:=Sqrt[x^2+y^2+z^2];

(*
Print["\[Alpha] = " <> ToString[\[Alpha],CForm]];
*)

(*
Print["\[Beta] = " <> ToString[\[Beta],CForm]];
*)

(*
ClearAll[i,j];
For[i=1,i<=3,i++,
  For[j=1,j\[LessEqual]3,j++,
    Print["\[Gamma]["<>ToString[i-1]<>"]["<>ToString[j-1]<>"] = "<>ToString[\[Gamma][[i,j]],CForm]];
  ]
]
*)

(*
ClearAll[i,j];
For[i=1,i<=3,i++,
  For[j=1,j\[LessEqual]3,j++,
    Print["i\[Gamma]["<>ToString[i-1]<>"]["<>ToString[j-1]<>"] = "<>ToString[i\[Gamma][[i,j]],CForm]];
  ]
]
*)

(*
ClearAll[i];
For[i=1,i<=3,i++,
Print["dbeta["<>ToString[i-1]<>"] = "<>ToString[gradi\[Beta]1]];
]
*)

(*
ClearAll[i];
For[i=1,i<=3,i++,
Print["dalpha["<>ToString[i-1]<>"] = "<>ToString[grad\[Alpha][[i]],CForm]];
]
*)

(*
ClearAll[i,j,k]
For[i=1,i\[LessEqual]3,i++,
  For[j=1,j\[LessEqual]3,j++,
    For[k=1,k\[LessEqual]3,k++,
    If[j\[LessEqual]k,Print["Gamma["<>ToString[i-1]<>"]["<>ToString[j-1]<>"]["<>ToString[k-1]<>"]= "<>ToString[\[CapitalGamma][i,j,k],CForm]]];
    ]
  ]
]
*)

ClearAll[r];



