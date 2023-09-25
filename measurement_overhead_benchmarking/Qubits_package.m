(* ::Package:: *)

(*Parameters*)


PI=N[\[Pi],12];
EPS=1.*^-12;


(*Pauli operators*)


\[Sigma]I={{1.,0.},{0.,1.}};
\[Sigma]X={{0.,1.},{1.,0.}};
\[Sigma]Y={{0.,-1.I},{1.I,0.}};
\[Sigma]Z={{1.,0.},{0.,-1.}};


funPX[i_]:=Module[{PX},(
	If[i==0,PX=\[Sigma]X,PX=\[Sigma]I];
	Do[(
		If[i==k,PX=KroneckerProduct[\[Sigma]X,PX],PX=KroneckerProduct[\[Sigma]I,PX]];
	),{k,1,Nq-1}];
	Return[PX]
)]

funPY[i_]:=Module[{PY},(
	If[i==0,PY=\[Sigma]Y,PY=\[Sigma]I];
	Do[(
		If[i==k,PY=KroneckerProduct[\[Sigma]Y,PY],PY=KroneckerProduct[\[Sigma]I,PY]];
	),{k,1,Nq-1}];
	Return[PY]
)]

funPZ[i_]:=Module[{PZ},(
	If[i==0,PZ=\[Sigma]Z,PZ=\[Sigma]I];
	Do[(
		If[i==k,PZ=KroneckerProduct[\[Sigma]Z,PZ],PZ=KroneckerProduct[\[Sigma]I,PZ]];
	),{k,1,Nq-1}];
	Return[PZ]
)]


funPO[ps_]:=Module[{q,p,PO},(
	q=0;
	p=StringTake[ps,{q+1}];
	If[StringMatchQ[p,"I"],PO=\[Sigma]I];
	If[StringMatchQ[p,"X"],PO=\[Sigma]X];
	If[StringMatchQ[p,"Y"],PO=\[Sigma]Y];
	If[StringMatchQ[p,"Z"],PO=\[Sigma]Z];
	Do[(
		p=StringTake[ps,{q+1}];
		If[StringMatchQ[p,"I"],PO=KroneckerProduct[\[Sigma]I,PO]];
		If[StringMatchQ[p,"X"],PO=KroneckerProduct[\[Sigma]X,PO]];
		If[StringMatchQ[p,"Y"],PO=KroneckerProduct[\[Sigma]Y,PO]];
		If[StringMatchQ[p,"Z"],PO=KroneckerProduct[\[Sigma]Z,PO]];
	),{q,1,Nq-1}];
	Return[PO]
)]


(*Hamiltonian*)


funHamiltonianQubit[Model_]:=Module[{Ham},(
	Ham=0.;
	Do[(
		Ham=Ham+Model[[i,1]]*funPO[Model[[i,2]]];
	),{i,1,Length[Model]}];
	Return[Ham]
)]


(*Spectrum*)


funSpectrum[Ham_]:=Module[{vals,vecs},(
	{vals,vecs}=Eigensystem[Ham];
	vals=Re[vals];
	{vals,vecs}=Transpose@SortBy[Transpose[{vals,vecs}],First];
	(*Print[Total[Total[Abs[Transpose[vecs].DiagonalMatrix[vals].Conjugate[vecs]-Ham]]]];*)
	(*Print[Total[Total[Abs[DiagonalMatrix[vals]-Conjugate[vecs].Ham.Transpose[vecs]]]]];*)
	Return[{vals,vecs}]
)]
