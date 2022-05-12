(* ::Package:: *)

(* ::Section:: *)
(*Spacetime metric (User)*)


(* The coordinate basis *)
ClearAll[coords];
coords = {t, x, y, z};

(* Auxiliary metric functions *)
ClearAll[lA, lB, \[ScriptCapitalH]A, \[ScriptCapitalH]B];
lA = {1, l1A[XA[t,x,y,z], YA[t,x,y,z], ZA[t,x,y,z]], l2A[XA[t,x,y,z], YA[t,x,y,z], ZA[t,x,y,z]], l3A[XA[t,x,y,z], YA[t,x,y,z], ZA[t,x,y,z]] };
lB = {1, l1B[XB[t,x,y,z], YB[t,x,y,z], ZB[t,x,y,z]], l2B[XB[t,x,y,z], YB[t,x,y,z], ZB[t,x,y,z]], l3B[XB[t,x,y,z], YB[t,x,y,z], ZB[t,x,y,z]] };
\[ScriptCapitalH]A = Table[2 * HA[XA[t,x,y,z], YA[t,x,y,z], ZA[t,x,y,z]] * lA[[a]] * lA[[b]], {a, 1, 4}, {b, 1,4}];
\[ScriptCapitalH]B = Table[2 * HB[XB[t,x,y,z], YB[t,x,y,z], ZB[t,x,y,z]] * lB[[a]] * lB[[b]], {a, 1, 4}, {b, 1,4}];

ClearAll[JA,JB];
JA = Table[D[{TA[t,x,y,z], XA[t,x,y,z], YA[t,x,y,z], ZA[t,x,y,z]}[[i]],coords[[j]]], {i, 1, 4}, {j, 1, 4}];
JB = Table[D[{TB[t,x,y,z], XB[t,x,y,z], YB[t,x,y,z], ZB[t,x,y,z]}[[i]],coords[[j]]], {i, 1, 4}, {j, 1, 4}];

ClearAll[eta];
eta = DiagonalMatrix[{-1, 1, 1, 1}];

ClearAll[T];
T[bhIdx_,t_,x_,y_,z_]={
{T[0,0,bhIdx,t,x,y,z],T[0,1,bhIdx,t,x,y,z],T[0,2,bhIdx,t,x,y,z],T[0,3,bhIdx,t,x,y,z]},
{T[0,1,bhIdx,t,x,y,z],T[1,1,bhIdx,t,x,y,z],T[1,2,bhIdx,t,x,y,z],T[1,3,bhIdx,t,x,y,z]},
{T[0,2,bhIdx,t,x,y,z],T[1,2,bhIdx,t,x,y,z],T[2,2,bhIdx,t,x,y,z],T[2,3,bhIdx,t,x,y,z]},
{T[0,3,bhIdx,t,x,y,z],T[1,3,bhIdx,t,x,y,z],T[2,3,bhIdx,t,x,y,z],T[3,3,bhIdx,t,x,y,z]}
};

ClearAll[ll4metric];
(*ll4metric = FullSimplify[
  Table[
    eta[[a, b]] + 
    Sum[JA[[ab,a]] * JA[[bb,b]] * \[ScriptCapitalH]A[[ab, bb]], {ab, 1, 4}, {bb, 1, 4}] +
    Sum[JB[[ab,a]] * JB[[bb,b]] * \[ScriptCapitalH]B[[ab, bb]], {ab, 1, 4}, {bb, 1, 4}],
    {a, 1, 4},{b, 1, 4}
  ]
];*)
ll4metric = eta + T[0,t,x,y,z]+T[1,t,x,y,z];


lapse


(* ::Section:: *)
(*ADM components*)


(* Assume that the coordinates are all real *)
$Assumptions = $Assumptions && And@@(Element[#,Reals]&/@coords);

(* The inverse metric. The user might want to set this directly for speed *)
ClearAll[uu4metric];
uu4metric = Inverse[ll4metric];

(* Spatial coordinates *)
ClearAll[spaceCoords];
spaceCoords = coords[[2;;4]];

ClearAll[llsmetric,uusmetric];
llsmetric = FullSimplify[ll4metric[[2;;4, 2;;4]]];
uusmetric = Inverse[llsmetric];

ClearAll[lapse];
lapse = Sqrt[-1/uu4metric[[1,1]]];

ClearAll[lshift,ushift];
lshift = ll4metric[[1, 2;;4]];
ushift = lapse^2 * uu4metric[[1, 2;;4]];

ClearAll[\[CapitalGamma],Dd];
\[CapitalGamma][i_,j_,k_] := FullSimplify[Sum[1/2 * uusmetric[[i,l]] * ( D[llsmetric[[l,j]],spaceCoords[[k]]] + D[llsmetric[[l,k]],spaceCoords[[j]]] - D[llsmetric[[j,k]],spaceCoords[[l]]] ),{l,1,3}]];
Dd[b_,c_,f_] := FullSimplify[D[f[[c]], spaceCoords[[b]]] - Sum[\[CapitalGamma][d,b,c]f[[d]],{d,1,3}]];

ClearAll[llextrinsic,ulextrinsic];
llextrinsic = FullSimplify[1/(2*lapse)*Table[-D[llsmetric[[i,j]],coords[[1]]] + Dd[i,j,lshift] + Dd[j,i,lshift], {i,1,3}, {j, 1,3}]];
ulextrinsic = FullSimplify[Table[Sum[uusmetric[[i,k]]*llextrinsic[[k,j]],{k,1,3}],{i,1,3},{j,1,3}]];


(* ::Section::Closed:: *)
(*ADM Derivatives*)


ClearAll[gradlapse];
gradlapse = Table[D[lapse,spaceCoords[[i]]], {i,1,3}];

ClearAll[gradushift];
gradushift = Table[D[ushift[[j]],spaceCoords[[i]]],{i,1,3},{j,1,3}];


(* ::Section:: *)
(*Substitution Rules*)


rules = {
HA[XA[t,x,y,z], YA[t,x,y,z], ZA[t,x,y,z]]->HA,
HB[XB[t,x,y,z], YB[t,x,y,z], ZB[t,x,y,z]]->HB,

l1A[XA[t,x,y,z], YA[t,x,y,z], ZA[t,x,y,z]]->l1A,
l2A[XA[t,x,y,z], YA[t,x,y,z], ZA[t,x,y,z]]->l2A,
l3A[XA[t,x,y,z], YA[t,x,y,z], ZA[t,x,y,z]]->l3A,

l1B[XB[t,x,y,z], YB[t,x,y,z], ZB[t,x,y,z]]->l1B,
l2B[XB[t,x,y,z], YB[t,x,y,z], ZB[t,x,y,z]]->l2B,
l3B[XB[t,x,y,z], YB[t,x,y,z], ZB[t,x,y,z]]->l3B,

D[HA[XA[t,x,y,z], YA[t,x,y,z], ZA[t,x,y,z]], XA[t,x,y,z]]->dHAdXA,
D[HA[XA[t,x,y,z], YA[t,x,y,z], ZA[t,x,y,z]], YA[t,x,y,z]]->dHAdYA,
D[HA[XA[t,x,y,z], YA[t,x,y,z], ZA[t,x,y,z]], ZA[t,x,y,z]]->dHAdZA,

D[HB[XB[t,x,y,z], YB[t,x,y,z], ZB[t,x,y,z]], XB[t,x,y,z]]->dHBdXB,
D[HB[XB[t,x,y,z], YB[t,x,y,z], ZB[t,x,y,z]], YB[t,x,y,z]]->dHBdYB,
D[HB[XB[t,x,y,z], YB[t,x,y,z], ZB[t,x,y,z]], ZB[t,x,y,z]]->dHBdZB,

D[l1A[XA[t,x,y,z], YA[t,x,y,z], ZA[t,x,y,z]], XA[t,x,y,z]]->dl1AdXA,
D[l1A[XA[t,x,y,z], YA[t,x,y,z], ZA[t,x,y,z]], YA[t,x,y,z]]->dl1AdYA,
D[l1A[XA[t,x,y,z], YA[t,x,y,z], ZA[t,x,y,z]], ZA[t,x,y,z]]->dl1AdZA,

D[l2A[XA[t,x,y,z], YA[t,x,y,z], ZA[t,x,y,z]], XA[t,x,y,z]]->dl2AdXA,
D[l2A[XA[t,x,y,z], YA[t,x,y,z], ZA[t,x,y,z]], YA[t,x,y,z]]->dl2AdYA,
D[l2A[XA[t,x,y,z], YA[t,x,y,z], ZA[t,x,y,z]], ZA[t,x,y,z]]->dl2AdZA,

D[l3A[XA[t,x,y,z], YA[t,x,y,z], ZA[t,x,y,z]], XA[t,x,y,z]]->dl3AdXA,
D[l3A[XA[t,x,y,z], YA[t,x,y,z], ZA[t,x,y,z]], YA[t,x,y,z]]->dl3AdYA,
D[l3A[XA[t,x,y,z], YA[t,x,y,z], ZA[t,x,y,z]], ZA[t,x,y,z]]->dl3AdZA,

D[l1B[XB[t,x,y,z], YB[t,x,y,z], ZB[t,x,y,z]], XB[t,x,y,z]]->dl1BdXB,
D[l1B[XB[t,x,y,z], YB[t,x,y,z], ZB[t,x,y,z]], YB[t,x,y,z]]->dl1BdYB,
D[l1B[XB[t,x,y,z], YB[t,x,y,z], ZB[t,x,y,z]], ZB[t,x,y,z]]->dl1BdZB,

D[l2B[XB[t,x,y,z], YB[t,x,y,z], ZB[t,x,y,z]], XB[t,x,y,z]]->dl2BdXB,
D[l2B[XB[t,x,y,z], YB[t,x,y,z], ZB[t,x,y,z]], YB[t,x,y,z]]->dl2BdYB,
D[l2B[XB[t,x,y,z], YB[t,x,y,z], ZB[t,x,y,z]], ZB[t,x,y,z]]->dl2BdZB,

D[l3B[XB[t,x,y,z], YB[t,x,y,z], ZB[t,x,y,z]], XB[t,x,y,z]]->dl3BdXB,
D[l3B[XB[t,x,y,z], YB[t,x,y,z], ZB[t,x,y,z]], YB[t,x,y,z]]->dl3BdYB,
D[l3B[XB[t,x,y,z], YB[t,x,y,z], ZB[t,x,y,z]], ZB[t,x,y,z]]->dl3BdZB,

D[TA[t,x,y,z],t]->jacA[0,0],
D[TA[t,x,y,z],x]->jacA[0,1],
D[TA[t,x,y,z],y]->jacA[0,2],
D[TA[t,x,y,z],z]->0(*jacA[0,3]*),

D[XA[t,x,y,z],t]->jacA[1,0],
D[XA[t,x,y,z],x]->jacA[1,1],
D[XA[t,x,y,z],y]->jacA[1,2],
D[XA[t,x,y,z],z]->0(*jacA[1,3]*),

D[YA[t,x,y,z],t]->jacA[2,0],
D[YA[t,x,y,z],x]->jacA[2,1],
D[YA[t,x,y,z],y]->jacA[2,2],
D[YA[t,x,y,z],z]->0(*jacA[2,3]*),

D[ZA[t,x,y,z],t]->0(*jacA[3,0]*),
D[ZA[t,x,y,z],x]->0(*jacA[3,1]*),
D[ZA[t,x,y,z],y]->0(*jacA[3,2]*),
D[ZA[t,x,y,z],z]->1(*jacA[3,3]*),

D[TB[t,x,y,z],t]->jacB[0,0],
D[TB[t,x,y,z],x]->jacB[0,1],
D[TB[t,x,y,z],y]->jacB[0,2],
D[TB[t,x,y,z],z]->0(*jacB[0,3]*),

D[XB[t,x,y,z],t]->jacB[1,0],
D[XB[t,x,y,z],x]->jacB[1,1],
D[XB[t,x,y,z],y]->jacB[1,2],
D[XB[t,x,y,z],z]->0(*jacB[1,3]*),

D[YB[t,x,y,z],t]->jacB[2,0],
D[YB[t,x,y,z],x]->jacB[2,1],
D[YB[t,x,y,z],y]->jacB[2,2],
D[YB[t,x,y,z],z]->0(*jacB[2,3]*),

D[ZB[t,x,y,z],t]->0(*jacB[3,0]*),
D[ZB[t,x,y,z],x]->0(*jacB[3,1]*),
D[ZB[t,x,y,z],y]->0(*jacB[3,2]*),
D[ZB[t,x,y,z],z]->1(*jacB[3,3]*),

D[D[TA[t,x,y,z],t],t]->djacA[0,0,0],
D[D[TA[t,x,y,z],t],x]->djacA[0,0,1],
D[D[TA[t,x,y,z],t],y]->djacA[0,0,2],
D[D[TA[t,x,y,z],t],z]->djacA[0,0,3],

D[D[TA[t,x,y,z],x],t]->djacA[0,1,0],
D[D[TA[t,x,y,z],x],x]->djacA[0,1,1],
D[D[TA[t,x,y,z],x],y]->djacA[0,1,2],
D[D[TA[t,x,y,z],x],z]->djacA[0,1,3],

D[D[TA[t,x,y,z],y],t]->djacA[0,2,0],
D[D[TA[t,x,y,z],y],x]->djacA[0,2,1],
D[D[TA[t,x,y,z],y],y]->djacA[0,2,2],
D[D[TA[t,x,y,z],y],z]->djacA[0,2,3],

D[D[TA[t,x,y,z],z],t]->djacA[0,3,0],
D[D[TA[t,x,y,z],z],x]->djacA[0,3,1],
D[D[TA[t,x,y,z],z],y]->djacA[0,3,2],
D[D[TA[t,x,y,z],z],z]->djacA[0,3,3],

D[D[XA[t,x,y,z],t],t]->djacA[1,0,0],
D[D[XA[t,x,y,z],t],x]->djacA[1,0,1],
D[D[XA[t,x,y,z],t],y]->djacA[1,0,2],
D[D[XA[t,x,y,z],t],z]->djacA[1,0,3],

D[D[XA[t,x,y,z],x],t]->djacA[1,1,0],
D[D[XA[t,x,y,z],x],x]->djacA[1,1,1],
D[D[XA[t,x,y,z],x],y]->djacA[1,1,2],
D[D[XA[t,x,y,z],x],z]->djacA[1,1,3],

D[D[XA[t,x,y,z],y],t]->djacA[1,2,0],
D[D[XA[t,x,y,z],y],x]->djacA[1,2,1],
D[D[XA[t,x,y,z],y],y]->djacA[1,2,2],
D[D[XA[t,x,y,z],y],z]->djacA[1,2,3],

D[D[XA[t,x,y,z],z],t]->djacA[1,3,0],
D[D[XA[t,x,y,z],z],x]->djacA[1,3,1],
D[D[XA[t,x,y,z],z],y]->djacA[1,3,2],
D[D[XA[t,x,y,z],z],z]->djacA[1,3,3],

D[D[YA[t,x,y,z],t],t]->djacA[2,0,0],
D[D[YA[t,x,y,z],t],x]->djacA[2,0,1],
D[D[YA[t,x,y,z],t],y]->djacA[2,0,2],
D[D[YA[t,x,y,z],t],z]->djacA[2,0,3],

D[D[YA[t,x,y,z],x],t]->djacA[2,1,0],
D[D[YA[t,x,y,z],x],x]->djacA[2,1,1],
D[D[YA[t,x,y,z],x],y]->djacA[2,1,2],
D[D[YA[t,x,y,z],x],z]->djacA[2,1,3],

D[D[YA[t,x,y,z],y],t]->djacA[2,2,0],
D[D[YA[t,x,y,z],y],x]->djacA[2,2,1],
D[D[YA[t,x,y,z],y],y]->djacA[2,2,2],
D[D[YA[t,x,y,z],y],z]->djacA[2,2,3],

D[D[YA[t,x,y,z],z],t]->djacA[2,3,0],
D[D[YA[t,x,y,z],z],x]->djacA[2,3,1],
D[D[YA[t,x,y,z],z],y]->djacA[2,3,2],
D[D[YA[t,x,y,z],z],z]->djacA[2,3,3],

D[D[ZA[t,x,y,z],t],t]->djacA[3,0,0],
D[D[ZA[t,x,y,z],t],x]->djacA[3,0,1],
D[D[ZA[t,x,y,z],t],y]->djacA[3,0,2],
D[D[ZA[t,x,y,z],t],z]->djacA[3,0,3],

D[D[ZA[t,x,y,z],x],t]->djacA[3,1,0],
D[D[ZA[t,x,y,z],x],x]->djacA[3,1,1],
D[D[ZA[t,x,y,z],x],y]->djacA[3,1,2],
D[D[ZA[t,x,y,z],x],z]->djacA[3,1,3],

D[D[ZA[t,x,y,z],y],t]->djacA[3,2,0],
D[D[ZA[t,x,y,z],y],x]->djacA[3,2,1],
D[D[ZA[t,x,y,z],y],y]->djacA[3,2,2],
D[D[ZA[t,x,y,z],y],z]->djacA[3,2,3],

D[D[ZA[t,x,y,z],z],t]->djacA[3,3,0],
D[D[ZA[t,x,y,z],z],x]->djacA[3,3,1],
D[D[ZA[t,x,y,z],z],y]->djacA[3,3,2],
D[D[ZA[t,x,y,z],z],z]->djacA[3,3,3],

D[D[TB[t,x,y,z],t],t]->djacB[0,0,0],
D[D[TB[t,x,y,z],t],x]->djacB[0,0,1],
D[D[TB[t,x,y,z],t],y]->djacB[0,0,2],
D[D[TB[t,x,y,z],t],z]->djacB[0,0,3],

D[D[TB[t,x,y,z],x],t]->djacB[0,1,0],
D[D[TB[t,x,y,z],x],x]->djacB[0,1,1],
D[D[TB[t,x,y,z],x],y]->djacB[0,1,2],
D[D[TB[t,x,y,z],x],z]->djacB[0,1,3],

D[D[TB[t,x,y,z],y],t]->djacB[0,2,0],
D[D[TB[t,x,y,z],y],x]->djacB[0,2,1],
D[D[TB[t,x,y,z],y],y]->djacB[0,2,2],
D[D[TB[t,x,y,z],y],z]->djacB[0,2,3],

D[D[TB[t,x,y,z],z],t]->djacB[0,3,0],
D[D[TB[t,x,y,z],z],x]->djacB[0,3,1],
D[D[TB[t,x,y,z],z],y]->djacB[0,3,2],
D[D[TB[t,x,y,z],z],z]->djacB[0,3,3],

D[D[XB[t,x,y,z],t],t]->djacB[1,0,0],
D[D[XB[t,x,y,z],t],x]->djacB[1,0,1],
D[D[XB[t,x,y,z],t],y]->djacB[1,0,2],
D[D[XB[t,x,y,z],t],z]->djacB[1,0,3],

D[D[XB[t,x,y,z],x],t]->djacB[1,1,0],
D[D[XB[t,x,y,z],x],x]->djacB[1,1,1],
D[D[XB[t,x,y,z],x],y]->djacB[1,1,2],
D[D[XB[t,x,y,z],x],z]->djacB[1,1,3],

D[D[XB[t,x,y,z],y],t]->djacB[1,2,0],
D[D[XB[t,x,y,z],y],x]->djacB[1,2,1],
D[D[XB[t,x,y,z],y],y]->djacB[1,2,2],
D[D[XB[t,x,y,z],y],z]->djacB[1,2,3],

D[D[XB[t,x,y,z],z],t]->djacB[1,3,0],
D[D[XB[t,x,y,z],z],x]->djacB[1,3,1],
D[D[XB[t,x,y,z],z],y]->djacB[1,3,2],
D[D[XB[t,x,y,z],z],z]->djacB[1,3,3],

D[D[YB[t,x,y,z],t],t]->djacB[2,0,0],
D[D[YB[t,x,y,z],t],x]->djacB[2,0,1],
D[D[YB[t,x,y,z],t],y]->djacB[2,0,2],
D[D[YB[t,x,y,z],t],z]->djacB[2,0,3],

D[D[YB[t,x,y,z],x],t]->djacB[2,1,0],
D[D[YB[t,x,y,z],x],x]->djacB[2,1,1],
D[D[YB[t,x,y,z],x],y]->djacB[2,1,2],
D[D[YB[t,x,y,z],x],z]->djacB[2,1,3],

D[D[YB[t,x,y,z],y],t]->djacB[2,2,0],
D[D[YB[t,x,y,z],y],x]->djacB[2,2,1],
D[D[YB[t,x,y,z],y],y]->djacB[2,2,2],
D[D[YB[t,x,y,z],y],z]->djacB[2,2,3],

D[D[YB[t,x,y,z],z],t]->djacB[2,3,0],
D[D[YB[t,x,y,z],z],x]->djacB[2,3,1],
D[D[YB[t,x,y,z],z],y]->djacB[2,3,2],
D[D[YB[t,x,y,z],z],z]->djacB[2,3,3],

D[D[ZB[t,x,y,z],t],t]->djacB[3,0,0],
D[D[ZB[t,x,y,z],t],x]->djacB[3,0,1],
D[D[ZB[t,x,y,z],t],y]->djacB[3,0,2],
D[D[ZB[t,x,y,z],t],z]->djacB[3,0,3],

D[D[ZB[t,x,y,z],x],t]->djacB[3,1,0],
D[D[ZB[t,x,y,z],x],x]->djacB[3,1,1],
D[D[ZB[t,x,y,z],x],y]->djacB[3,1,2],
D[D[ZB[t,x,y,z],x],z]->djacB[3,1,3],

D[D[ZB[t,x,y,z],y],t]->djacB[3,2,0],
D[D[ZB[t,x,y,z],y],x]->djacB[3,2,1],
D[D[ZB[t,x,y,z],y],y]->djacB[3,2,2],
D[D[ZB[t,x,y,z],y],z]->djacB[3,2,3],

D[D[ZB[t,x,y,z],z],t]->djacB[3,3,0],
D[D[ZB[t,x,y,z],z],x]->djacB[3,3,1],
D[D[ZB[t,x,y,z],z],y]->djacB[3,3,2],
D[D[ZB[t,x,y,z],z],z]->djacB[3,3,3]
};


(* ::Section:: *)
(*Expression output*)


Export["lapse.dat",ToString[lapse//.rules,CForm]<>";"];

Export["lshift_0.dat",ToString[lshift[[1]]//.rules,CForm]<>";"];
Export["lshift_1.dat",ToString[lshift[[2]]//.rules,CForm]<>";"];
Export["lshift_2.dat",ToString[lshift[[3]]//.rules,CForm]<>";"];

Export["ushift_0.dat",ToString[ushift[[1]]//.rules,CForm]<>";"];
Export["ushift_1.dat",ToString[ushift[[2]]//.rules,CForm]<>";"];
Export["ushift_2.dat",ToString[ushift[[3]]//.rules,CForm]<>";"];

Export["llsmetric_00.dat",ToString[llsmetric[[1,1]]//.rules,CForm]<>";"];
Export["llsmetric_01.dat",ToString[llsmetric[[1,2]]//.rules,CForm]<>";"];
Export["llsmetric_02.dat",ToString[llsmetric[[1,3]]//.rules,CForm]<>";"];
Export["llsmetric_11.dat",ToString[llsmetric[[2,2]]//.rules,CForm]<>";"];
Export["llsmetric_12.dat",ToString[llsmetric[[2,3]]//.rules,CForm]<>";"];
Export["llsmetric_22.dat",ToString[llsmetric[[3,3]]//.rules,CForm]<>";"];

Export["uusmetric_00.dat",ToString[uusmetric[[1,1]]//.rules,CForm]<>";"];
Export["uusmetric_01.dat",ToString[uusmetric[[1,2]]//.rules,CForm]<>";"];
Export["uusmetric_02.dat",ToString[uusmetric[[1,3]]//.rules,CForm]<>";"];
Export["uusmetric_11.dat",ToString[uusmetric[[2,2]]//.rules,CForm]<>";"];
Export["uusmetric_12.dat",ToString[uusmetric[[2,3]]//.rules,CForm]<>";"];
Export["uusmetric_22.dat",ToString[uusmetric[[3,3]]//.rules,CForm]<>";"];

Export["llextrinsic_00.dat",ToString[llextrinsic[[1,1]]//.rules,CForm]<>";"];
Export["llextrinsic_01.dat",ToString[llextrinsic[[1,2]]//.rules,CForm]<>";"];
Export["llextrinsic_02.dat",ToString[llextrinsic[[1,3]]//.rules,CForm]<>";"];
Export["llextrinsic_11.dat",ToString[llextrinsic[[2,2]]//.rules,CForm]<>";"];
Export["llextrinsic_12.dat",ToString[llextrinsic[[2,3]]//.rules,CForm]<>";"];
Export["llextrinsic_22.dat",ToString[llextrinsic[[3,3]]//.rules,CForm]<>";"];

Export["ulextrinsic_00.dat",ToString[ulextrinsic[[1,1]]//.rules,CForm]<>";"];
Export["ulextrinsic_01.dat",ToString[ulextrinsic[[1,2]]//.rules,CForm]<>";"];
Export["ulextrinsic_02.dat",ToString[ulextrinsic[[1,3]]//.rules,CForm]<>";"];
Export["ulextrinsic_11.dat",ToString[ulextrinsic[[2,2]]//.rules,CForm]<>";"];
Export["ulextrinsic_12.dat",ToString[ulextrinsic[[2,3]]//.rules,CForm]<>";"];
Export["ulextrinsic_22.dat",ToString[ulextrinsic[[3,3]]//.rules,CForm]<>";"];

Export["gradlapse_0.dat",ToString[gradlapse[[1]]//.rules,CForm]<>";"];
Export["gradlapse_1.dat",ToString[gradlapse[[2]]//.rules,CForm]<>";"];
Export["gradlapse_2.dat",ToString[gradlapse[[3]]//.rules,CForm]<>";"];

Export["gradushift_00.dat",ToString[gradushift[[1,1]]//.rules,CForm]<>";"];
Export["gradushift_01.dat",ToString[gradushift[[1,2]]//.rules,CForm]<>";"];
Export["gradushift_02.dat",ToString[gradushift[[1,3]]//.rules,CForm]<>";"];
Export["gradushift_11.dat",ToString[gradushift[[2,2]]//.rules,CForm]<>";"];
Export["gradushift_12.dat",ToString[gradushift[[2,3]]//.rules,CForm]<>";"];
Export["gradushift_22.dat",ToString[gradushift[[3,3]]//.rules,CForm]<>";"];
