(* ::Package:: *)

BeginPackage["VQE`"];

CircleTimes::usage = "Standard KroneckerProduct function. It can deal \
with indices. Assumes ordered matrices, from low index to high. \
\!\(\*SubscriptBox[\(A\), \(2\)]\)\[CircleTimes]\!\(\*SubscriptBox[\(B\), \(3\)]\)\[CircleTimes]\!\(\*SubscriptBox[\(C\), \(5\)]\) = KroneckerProduct[A,B,C].";
CircleDot::usage = "Equivalent to CircleTimes function, but also pads \
the resulting matrix, to account for \!\(\*SubscriptBox[\(N\), \(qubits\)]\). \
\!\(\*SubscriptBox[\(A\), \(2\)]\)\[CircleDot]\!\(\*SubscriptBox[\(B\), \(3\)]\)\[CircleDot]\!\(\*SubscriptBox[\(C\), \(5\)]\) = \!\(\*SubscriptBox[\(Id\), \(1\)]\)\[CircleTimes]\!\(\*SubscriptBox[\(A\), \(2\)]\)\[CircleTimes]\!\(\*SubscriptBox[\(B\), \(3\)]\)\[CircleTimes]\!\(\*SubscriptBox[\(Id\), \(4\)]\)\[CircleTimes]\!\(\*SubscriptBox[\(C\), \(5\)]\)\[CircleTimes]\!\(\*SubscriptBox[\(Id\), \(6\)]\) (if \!\(\*SubscriptBox[\(N\), \(qubits\)]\)=6).";

Nqubits::usage = "\!\(\*SubscriptBox[\(N\), \(qubits\)]\) in the system of study (default is 4).";
X::usage = "Pauli \!\(\*SuperscriptBox[\(\[Sigma]\), \(x\)]\) matrix.";
Y::usage = "Pauli \!\(\*SuperscriptBox[\(\[Sigma]\), \(y\)]\) matrix.";
Z::usage = "Pauli \!\(\*SuperscriptBox[\(\[Sigma]\), \(z\)]\) matrix.";
Id::usage = "2x2 Identity matrix.";

Hadg::usage = "Hadg[n] represents the Hadamard gate applied to qubit n \
(\!\(\*SubscriptBox[\(N\), \(qubits\)]\)\!\(\*SubscriptBox[\(xN\), \(qubits\)]\) matrix).";
CNOT::usage = "CNOT[c,t] represents the CNOT gate with control qubit c \
and target qubit t (\!\(\*SubscriptBox[\(N\), \(qubits\)]\)\!\(\*SubscriptBox[\(xN\), \(qubits\)]\) matrix).";
SWAP::usage = "CNOT[c,t] represents the CNOT gate with control qubit c \
and target qubit t (\!\(\*SubscriptBox[\(N\), \(qubits\)]\)\!\(\*SubscriptBox[\(xN\), \(qubits\)]\) matrix).";
Idg::usage = "Idg[n] represents the Identity gate of size n x n.";
Rot::usage = "Rot[M,\[Theta]] returns the \!\(\*SuperscriptBox[\(e\), \(\(-\*FractionBox[\(\[ImaginaryI]\), \(2\)]\) \(\[Theta]\)\(\\\ \)\(M\)\(\\\ \)\)]\)matrix.";
Pad::usage = "Pad[\!\(\*SubscriptBox[\(M\), \(n\)]\)], where M is a 2x2 matrix, returns \
\!\(\*SubscriptBox[\(Id\), \(1\)]\)\[CircleTimes]...\[CircleTimes]\!\(\*SubscriptBox[\(M\), \(n\)]\)\[CircleTimes]...\[CircleTimes]\!\(\*SubscriptBox[\(Id\), \(Nqubits\)]\)";

HtoM::usage = "Converts an expression written in second quantization \
to its equivalent matrix (using the Jordan-Wigner transformation).";
HtoBKM::usage = "Converts an expression written in second quantization \
to its equivalent matrix (using the Bravyi-Kitaev transformation).";
HtoS::usage = "Converts an expression written in second quantization \
to its equivalent simplified expression with Pauli matrices \
(using the Jordan-Wigner transformation).";
HtoP::usage = "Converts an expression written in second quantization \
to its equivalent Mathematica-usable expression with Pauli matrices \
(using the Jordan-Wigner transformation).";
PtoM::usage = "Converts a Mathematica-usable expression with Pauli \
matrices to its equivalent matrix. It is a placeholder function, equivalent to pressing \
Enter manually.";
PtoS::usage = "Converts a Mathematica-usable expression with Pauli \
matrices to its equivalent simplified expression with Pauli matrices.";
StoM::usage = "Converts a simplified expression with Pauli \
matrices to its equivalent matrix.";
StoP::usage = "Converts a simplified expression with Pauli \
matrices to its equivalent Mathematica-usable expression with Pauli \
matrices.";
MtoP::usage = "Converts a matrix to its equivalent Mathematica-usable \
expression with Pauli matrices.";
MtoS::usage = "Converts a matrix to its equivalent simplified \
expression with Pauli matrices.";

IniState::usage = "IniState[abcd] (abcd as a string) returns ket \
|abcd> as a column vector, with a being the lowest index qubit, and \
d the highest.";
FromVectorToState::usage = "FromVectorToState[x], with x being a \
column or line vector, converts x to its corresponding state in \
Dirac notation.";
Hypercube::usage = "Represents geometric region equivalent to \
Mathematica's Cuboid function, but easier to input for higher \
dimensions. Hypercube[{a,b},n] gives the cuboid defined from \
point {a,...,a} to point {b,...,b}, in n dimensions."


Begin["`Private`"];

Nqubits=4;

(*Operators*)

CircleTimes[x__]:=Module[{xa,xb},
xa=Replace[{x}//.{Subscript[y_, n_]-> y},y_/;Length@Dimensions[y]<= 1->{{y}},1];
xb=KroneckerProduct@@xa;
xb
];

CircleDot[x__]:=Module[{subs=Cases[{x},Subscript[a_, n_]-> {a,n}],mults},
mults=Table[If[MemberQ[subs[[All,2]], i],Cases[subs,{a_,i}-> a][[1]],Id],{i,1,Nqubits}];
KroneckerProduct@@mults
];

CirclePlus[x__]:=Module[{subs=Cases[{x},Subscript[a_, n_]-> {a,n}],mults,subsAfterList},
subsAfterList=Flatten[Table[If[IntegerQ[subs[[i,2]]],{subs[[i]]},Partition[Reverse[Riffle[subs[[i,2]],{subs[[i,1]]},{2,-1,2}]],2]],{i,1,Length@subs}],1];
mults=Table[If[MemberQ[subsAfterList[[All,2]], i],Cases[subsAfterList,{a_,i}-> a][[1]],Id],{i,1,Nqubits}];
If[Nqubits==1,Flatten[mults,1],KroneckerProduct@@mults]
];


(*Matrices*)
X={{0,1},{1,0}};
Y={{0,-I},{I,0}};
Z={{1,0},{0,-1}};
Id={{1,0},{0,1}};
\[Sigma]v[1]=Id;
\[Sigma]v[2]=X;
\[Sigma]v[3]=Y;
\[Sigma]v[4]=Z;
\[Sigma]vd[1]:=Id;
\[Sigma]vd[2]:=X;
\[Sigma]vd[3]:=Y;
\[Sigma]vd[4]:=Z;
\[Sigma]s[1]="Id";
\[Sigma]s[2]="X";
\[Sigma]s[3]="Y";
\[Sigma]s[4]="Z";

Had={{1,1},{1,-1}}/Sqrt[2];
CNOT2to1={{1,0,0,0},{0,1,0,0},{0,0,0,1},{0,0,1,0}};
CNOT0=(Had\[CircleTimes]Had).(CNOT2to1).(Had\[CircleTimes]Had);
Idg[n_]:=IdentityMatrix[2^n];
SWAPdef=CNOT0.CNOT2to1.CNOT0;
SWAPclose[a_]:=outerPad[SWAPdef,a-2,Nqubits-a];

SWAP[a_,b_]:=Module[{m,l,sl,in=Idg[Nqubits]},
If[a==b,in,
m=Max[a,b];
l=Min[a,b];
sl=If[l+1<m,Dot@@SWAPclose/@Range[l+1,m],in];
res=sl\[ConjugateTranspose].SWAPclose[l+1].sl]
];

CNOT[a_,b_]:=SWAP[a,2].SWAP[If[a==1&&b==2,1,b],1].rightPad[CNOT0,Nqubits-2].SWAP[1,If[a==1&&b==2,1,b]].SWAP[2,a];

Hadg[n_]:=Pad[Subscript[Had, n]];

Rot[m_,\[Theta]_]:=MatrixExp[-I \[Theta]/2 m];

Dim[m_]:=Log[2,Dimensions[m,1][[1]]];

leftPad[m_,n_]:=Idg[n]\[CircleTimes]m;
rightPad[m_,n_]:=m\[CircleTimes]Idg[n];
outerPad[m_,nl_,nr_]:=Idg[nl]\[CircleTimes]m\[CircleTimes]Idg[nr];
Pad[m_]:=(x=Cases[{m},Subscript[a_, n_]-> {a,n}][[1]];outerPad[x[[1]],x[[2]]-1,Nqubits-x[[2]]]);


(*Transformations*)
JW[a_]:=Module[{num=Cases[{a},Subscript[b_, n_]-> n][[1]],res},
res=KroneckerProduct@@Join[If[num==1,{{1}},ConstantArray[Z,num-1]],{((X+I Y)/2),Idg[Nqubits-num]}];
res
];

quadOne[d_]:=PadLeft[{ConstantArray[1,d]},d,{ConstantArray[0,d]}];
quadOnlyOne[d_]:=PadLeft[{PadLeft[{1},d,0]},d,{ConstantArray[0,d]}];
parity[no_]:=Table[If[i>j,1,0],{i,1,no},{j,1,no}];
beta[no_]:=Module[{no2=2^Ceiling[Log2[no]],finalBeta},
finalBeta=If[no2==1,{{1}},
({
 {beta[no2/2], 0},
 {quadOne[no2/2], beta[no2/2]}
})//ArrayFlatten];
finalBeta[[no2-no+1;;,no2-no+1;;]]
];
inverseBeta[n_]:=Mod[ModularInverse[Det[beta[n]],2] Det[beta[n]]Inverse[beta[n]],2];
inverseBetaAlt[no_]:=Module[{no2=2^Ceiling[Log2[no]],finalBeta},
finalBeta=If[no2==1,{{1}},
({
 {inverseBetaAlt[no2/2], 0},
 {quadOnlyOne[no2/2], inverseBetaAlt[no2/2]}
})//ArrayFlatten];
finalBeta[[no2-no+1;;,no2-no+1;;]]
];
inverseBetaUsed=inverseBetaAlt;

paritySet[n_]:=Module[{rowN,res},
rowN = Mod[parity[Nqubits].inverseBetaUsed[Nqubits],2][[n]];
res=Select[Position[rowN,_?(#==1&)]//Flatten,#<n& ];
res
];
updateSet[n_]:=Module[{colN,res},
colN = (beta[Nqubits]\[Transpose])[[n]];
res=Select[Position[colN,_?(#==1&)]//Flatten,#>n& ];
res
];
flipSet[n_]:=Module[{rowN,res},
rowN = inverseBetaUsed[Nqubits][[n]];
res=Select[Position[rowN,_?(#==1&)]//Flatten,#<n& ];
res
];
finalSet[n_]:=If[OddQ[n],paritySet[n],Complement[paritySet[n],flipSet[n]]];(*We use OddQ instead of EvenQ since this notation starts at 1 and not 0*)

BK[a_]:=Module[{num=Cases[{a},Subscript[b_, n_]->n][[1]],res},
res=1/2 (Subscript[X, updateSet[num]]\[CirclePlus]Subscript[X, num]\[CirclePlus]Subscript[Z, paritySet[num]]+ I*(Subscript[X, updateSet[num]]\[CirclePlus]Subscript[Y, num]\[CirclePlus]Subscript[Z, finalSet[num]]));
res
];


(*Conversions*)
HtoBKM[x__]:=ReplaceAll[x,Subscript[Global`a, n_]:>BK[Subscript[Global`a, n]]];
HtoM[x__]:=ReplaceAll[x,Subscript[Global`a, n_]:>JW[Subscript[Global`a, n]]];

\[Sigma]Decomposition[m_/;IntegerQ[Log[2,Length[m]]]] :=
  Module[{tiefe = Log[2,Length[m]],indextable,indexlist},
    indextable = Table[{{1},{2},{3},{4}},{Log[2,Length[m]]}];
    indexlist = Outer[Join,Sequence @@indextable,1];
    Map[\[Sigma]Extract[m,#]&, indexlist,{tiefe}]];
    
\[Sigma]Extract[m_,indices_]:=
  TrProd[m,CircleTimes[Sequence @@ \[Sigma]v /@ indices]]/Length[m];

TrProd[m1_,m2_]:=
	Sum[m1[[i,j]]*m2[[j,i]],{i,1,Length[m1]},{j,1,Length[m1]}];

MtoP[x_]:=Module[{tableSize},
tableSize=Module[{ts=ConstantArray[0,ConstantArray[4,Nqubits]],posTable,subsfun},
posTable=Position[ts,0];
subsfun[xo_]:=Subscript[\[Sigma]v[xo[[#]]],#]&/@Range[Nqubits];
ReplacePart[ts,# ->Inactive[CircleTimes][Sequence @@ subsfun[#]]&/@posTable]//.{{{1,0},{0,1}}:>Defer@Id,{{0,1},{1,0}}:> Defer@X,{{0,-I},{I,0}}:>Defer@ Y,{{1,0},{0,-1}}:> Defer@Z}
];
ToExpression[StringReplace[ToString@InputForm[(Simplify[\[Sigma]Decomposition[x]] tableSize)//Flatten//Total//FullSimplify],"Inactive[CircleTimes]"->"CircleTimes"],StandardForm,Defer]//FullSimplify
];

MtoS[x_]:=Module[{expr},
expr=MtoP[x];
expr=PtoS[expr]
];

PtoSold[x_]:=Module[{expr},
expr=x//.{HoldPattern[q__\[CircleTimes]w__]:>Global`qwak q Global`qwak w Global`qwak};Print[expr];
expr=expr//.{HoldPattern[q__ Subscript[Id, n_]w__]:>q w};Print[expr];
expr=ToString[expr,InputForm];Print[expr];
expr=StringReplace[StringReplace[expr,{"qwak*"-> "","*qwak"-> ""}],"Defer[qwak]"-> "1"];Print[expr];
expr=ToExpression[expr]
];
SetAttributes[PtoS,HoldAll];
PtoS[x_]:=Module[{expr,test,xo},
test=ToString[x,InputForm];
If[StringFreeQ[test,"Defer"],expr=StringReplace[ToString[HoldForm[x],InputForm],
{"Id":>"Defer[Id]","X":>"Defer[X]","Y":>"Defer[Y]","Z":>"Defer[Z]","HoldForm"->"Defer"}],
expr=test];
expr=StringReplace[expr," \[CircleTimes] "->"*"];
expr=StringReplace[StringReplace[expr, {"*Subscript[Defer[Id], "~~NumberString~~"]"->"",
"Subscript[Defer[Id], "~~NumberString~~"]*"->"","Subscript[Defer[Id], "~~NumberString~~"]"->""}],"*()"->""];
expr=StringDrop[StringDrop[expr,6],-1];
expr=ToExpression[expr,StandardForm]
];

SetAttributes[StoMA,HoldAll];
SetAttributes[StoMB,HoldAll];
SetAttributes[StoM,HoldAll];
StoMA[x_]:=Module[{xi,xa},
xi=ToString@InputForm@If[FreeQ[x//HoldForm//InputForm,Subscript],x,HoldForm@x];
xa=StringReplace[Nest[StringReplace[#,"Subscript["~~w_~~", "~~n_~~"]*Subscript["~~q_~~", "~~m_~~"]":>"Subscript["~~w~~", "~~n~~"]\[CircleDot]Subscript["~~q~~", "~~m~~"]"]&,xi,Nqubits-1],"Defer"-> "Identity"];
xa=If[StringFreeQ[xa,"HoldForm"],"HoldForm["<>xa<>"]",xa];
xa=ToExpression[xa,InputForm]
];
StoMB[x_]:=Module[{xa,expr},
xa=Inactivate[Evaluate@x,CircleDot];
xa=ExpandAll[ReleaseHold@xa];
expr=Table[Take[xa,{i}],{i,1,Length[xa]}];
expr=Table[If[FreeQ[expr[[i]]//FullForm,CircleDot],If[Cases[expr[[i]],Subscript[X, _]|Subscript[Y, _]|Subscript[Z, _]]=={},ReplaceAll[expr[[i]],q__->q Idg[Nqubits]],ReplaceAll[expr[[i]],{Subscript[X, n_]:>Pad[Subscript[X, n]],Subscript[Y, n_]:> Pad[Subscript[Y, n]],Subscript[Z, n_]:>Pad[Subscript[Z, n]] }]], ReplaceAll[expr[[i]],{q_\[CircleTimes]w_->q\[CircleDot]w}]],{i,1,Length[expr]}];
Total[expr//Activate]
];
StoM[x_]:=Module[{xa},
xa=StoMA[x];
xa=StoMB[xa]
];

StoP[x_]:=Module[{expr},
expr=StoM[x];
expr=MtoP[expr]
];

PtoM[x_]:=x//ToBoxes//ToExpression//FullSimplify;
HtoP[x_]:=x//HtoM//MtoP;
HtoS[x_]:=x//HtoM//MtoS;


(*Helpful functions*)
IniState[x_]:=Module[{sx=ToString[x],dec,vec},
dec=FromDigits[sx,2];
vec=ConstantArray[0,2^StringLength[sx]];
vec[[dec+1]]=1;
{vec}//Transpose
];

FromVectorToState[x_]:=Module[{xo,res},
If[Dimensions[x]=={2^Nqubits},xo=x,
If[Dimensions[x]=={1,2^Nqubits},xo=x[[1]],
If[Dimensions[x]=={2^Nqubits,1},xo=(x\[Transpose])[[1]],xo=IniState[StringPadRight["1",Nqubits,"0"]]]]];
res=Total@Table[xo[[i]]*("|"~~StringPadLeft[StringJoin[ToString/@IntegerDigits[i-1,2]],Nqubits,"0"]~~">"),{i,1,2^Nqubits}];
res//Chop
];

Hypercube[{a_,b_},n_]:=Cuboid[ConstantArray[a,n],ConstantArray[b,n]];

End[];

EndPackage[];
