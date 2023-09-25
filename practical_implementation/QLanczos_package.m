(* ::Package:: *)

(*Pauli operators and spins*)


\[Sigma]I={{1,0},{0,1}};
\[Sigma]X={{0,1},{1,0}};
\[Sigma]Y={{0,-I},{I,0}};
\[Sigma]Z={{1,0},{0,-1}};
spinUp={{1},{0}};
spinDown={{0},{1}};


(*Model and Hamiltonian*)


HeisenbergHam:=Module[{Ham,hTempX,hTempY,hTempZ,numQubit},(*J=1*)
numQubit=10;
Ham=ConstantArray[0,{2^numQubit,2^numQubit}];
Do[
   Do[
      If[(i==1 && j==2)||(i==2 && j==3)||(i==3 && j==4)||(i==4 && j==5)||(i==5 && j==6)||
         (i==6 && j==7)||(i==7 && j==8)||(i==8 && j==9)||(i==9 && j==10),(*lattice, i<j*)
         hTempX={{1}};hTempY={{1}};hTempZ={{1}};
         Do[
            If[k==i || k==j,
               hTempX=KroneckerProduct[hTempX,\[Sigma]X];hTempY=KroneckerProduct[hTempY,\[Sigma]Y];hTempZ=KroneckerProduct[hTempZ,\[Sigma]Z];
              ,hTempX=KroneckerProduct[hTempX,\[Sigma]I];hTempY=KroneckerProduct[hTempY,\[Sigma]I];hTempZ=KroneckerProduct[hTempZ,\[Sigma]I];
               ];
         ,{k,1,numQubit}];
         Ham=Ham+hTempX+hTempY+hTempZ;
      ];
   ,{j,1,numQubit}];
,{i,1,numQubit}];
Ham=N[Ham]]


HubbardHam[U_]:=Module[{Ham,J,hTempX,hTempY,hTempZ,numSite,numQubit},
J=1;
numSite=5;
numQubit=2*numSite;
Ham=ConstantArray[0,{2^numQubit,2^numQubit}];
(*interaction*)
Do[
hTempZ={{1}};
   Do[
      If[k==i||k==i+numSite,hTempZ=KroneckerProduct[hTempZ,\[Sigma]Z],hTempZ=KroneckerProduct[hTempZ,\[Sigma]I]];
   ,{k,1,numQubit}];
Ham=Ham+U/4 hTempZ;
,{i,1,numSite}];
(*hopping*)
Do[
   Do[
      If[(i==1 && j==2)||(i==1 && j==3)||(i==2 && j==3)||(i==2 && j==4)||(i==3 && j==4)||(i==3 && j==5)||(i==4 && j==5),(*lattice, i<j*)
         (*hopping_spinup*)
         hTempX={{1}};hTempY={{1}};
         If[i!=1,
            Do[
               hTempX=KroneckerProduct[hTempX,\[Sigma]I];hTempY=KroneckerProduct[hTempY,\[Sigma]I];
            ,{k,1,i-1}];
         ];
         Do[
            If[k==i||k==j,
                hTempX=KroneckerProduct[hTempX,\[Sigma]X];hTempY=KroneckerProduct[hTempY,\[Sigma]Y];
               ,hTempX=KroneckerProduct[hTempX,\[Sigma]Z];hTempY=KroneckerProduct[hTempY,\[Sigma]Z];];
         ,{k,i,j}];
         Do[
            hTempX=KroneckerProduct[hTempX,\[Sigma]I];hTempY=KroneckerProduct[hTempY,\[Sigma]I];
         ,{k,j+1,numQubit}];
         Ham=Ham-J/2 hTempX-J/2 hTempY;
         (*hopping_spindown*)
         hTempX={{1}};hTempY={{1}};
         Do[
            hTempX=KroneckerProduct[hTempX,\[Sigma]I];hTempY=KroneckerProduct[hTempY,\[Sigma]I];
         ,{k,1,i+numSite-1}];
         Do[
            If[k==i+numSite||k==j+numSite,
                hTempX=KroneckerProduct[hTempX,\[Sigma]X];hTempY=KroneckerProduct[hTempY,\[Sigma]Y];
               ,hTempX=KroneckerProduct[hTempX,\[Sigma]Z];hTempY=KroneckerProduct[hTempY,\[Sigma]Z];];
         ,{k,i+numSite,j+numSite}];
         If[j!=numSite,
            Do[
               hTempX=KroneckerProduct[hTempX,\[Sigma]I];hTempY=KroneckerProduct[hTempY,\[Sigma]I];
            ,{k,j+numSite+1,numQubit}];
         ];
         Ham=Ham-J/2 hTempX-J/2 hTempY;
      ];
   ,{j,1,numSite}];
,{i,1,numSite}];
Ham=N[Ham]]


(*Reference state*)


\[CurlyPhi]Heisenberg=(1./Sqrt[2])^5 (KroneckerProduct[
KroneckerProduct[spinUp,spinDown]-KroneckerProduct[spinDown,spinUp],
KroneckerProduct[spinUp,spinDown]-KroneckerProduct[spinDown,spinUp],
KroneckerProduct[spinUp,spinDown]-KroneckerProduct[spinDown,spinUp],
KroneckerProduct[spinUp,spinDown]-KroneckerProduct[spinDown,spinUp],
KroneckerProduct[spinUp,spinDown]-KroneckerProduct[spinDown,spinUp]]);


\[CurlyPhi]Hubbard:=Module[{Ham,vals,vecs,\[Psi]},
Ham=HubbardHam[0];
{vals,vecs}=funSpectrum[Ham];
\[Psi]=Transpose[{vecs[[1]]}];
\[Psi]]


(*Spectrum*)


funSpectrum[Ham_]:=Module[{vals,vecs},
{vals,vecs}=Eigensystem[Ham];
vals=Re[vals];
{vals,vecs}=Transpose@SortBy[Transpose[{vals,vecs}],First];
(*Total[Total[Abs[Transpose[vecs].DiagonalMatrix[vals].Conjugate[vecs]-Ham]]]*)
{vals,vecs}]


(*Subspace diagonalization*)


funSubDiag[Hmat_,Smat_]:=Module[{Svals,Svecs,V,Heff,vals,vecs,EK,cn},
{Svals,Svecs}=funSpectrum[Smat];
cn=Max[Abs[Svals]]/Min[Abs[Svals]];
V=Transpose[Svecs] . DiagonalMatrix[1./Sqrt[Svals]];
Heff=ConjugateTranspose[V] . Hmat . V; 
{vals,vecs}=funSpectrum[Heff];
EK=vals[[1]];
{EK,cn}]


(*(a^\dag S a)/(a^\dag a)*)


funaSa[Hmat_,Smat_]:=Module[{Svals,Svecs,V,Heff,vals,vecs,EK,cn,a,aSa},
{Svals,Svecs}=funSpectrum[Smat];
cn=Max[Abs[Svals]]/Min[Abs[Svals]];
V=Transpose[Svecs] . DiagonalMatrix[1./Sqrt[Svals]];
Heff=ConjugateTranspose[V] . Hmat . V; 
{vals,vecs}=funSpectrum[Heff];
EK=vals[[1]];
a=V . Transpose[vecs[[1]]];
aSa=Re[ConjugateTranspose[a] . Smat . a]/Re[ConjugateTranspose[a] . a];
{EK,cn,aSa}]


(*Cost for Gaussian power*)


funLORfactor[htot_,\[Tau]_,NT_,u_]:=Module[{t,factor},(
	t=\[Tau]*u/NT;
	factor=(Sqrt[1.+htot^2*t^2]+Exp[htot*t]-(1+htot*t))^NT/Exp[Exp[1.]*htot^2*t^2/2.*NT];
	Return[factor]
)]


funIntegral[htot_,\[Tau]_,d_]:=Module[{NT,\[Chi],costList,factor,cost},(
	NT=Ceiling[4.*Exp[1.]*htot^2*\[Tau]^2];
	\[Chi]=Exp[1.]*htot^2*\[Tau]^2/(2.*NT);
	costList={};
	Do[(
		If[n==0,factor=1.,factor=(n/Exp[1.])^(-(n/2))];
		cost=NIntegrate[2(2^(-(n/2)) Abs[HermiteH[n,u/Sqrt[2]]] 1/Sqrt[2\[Pi]] Exp[-(1/2-\[Chi]) u^2])*factor*funLORfactor[htot,\[Tau],NT,u],{u,0.,\[Infinity]},Method->"QuasiMonteCarlo"];
		AppendTo[costList,cost]
	),{n,0,d-1}];
	Return[costList]
)]


funCost[htot_,\[Tau]_,d_]:=Module[{costList},(
	costList=funIntegral[htot,\[Tau],d];
	Do[(
		costList[[k]]=costList[[k]]*(k-1)^((k-1)/2)/Exp[(k-1)/2]/\[Tau]^(k-1);
	),{k,2,d}];
	Return[costList]
)]


(*funCostGP[\[Tau]_,n_]:=Module[{\[Chi],A,cost},
\[Chi]=0.125;(*\[Chi]=1/8 for n=0 and \[Chi]=1/16 for n>0*)
A=1.537931098297999;
If[n==0,
   cost=Sqrt[1./(1.-2.*\[Chi])];,
   cost=A*n^(n/2)/Exp[n/2]/\[Tau]^n;
   ];
cost]*)


(*Hmat and Smat for different cases in Table I*)


(*1. Power*)


funMatP[\[CapitalLambda]_,E0_,d_,prob\[CurlyPhi]_]:=Module[{Hmat,Smat},
Hmat=ConstantArray[0,{d,d}];
Smat=ConstantArray[0,{d,d}];
Do[
   Do[
      If[k+q-2==0,
         Hmat[[k,q]]=Total[prob\[CurlyPhi]*\[CapitalLambda]];
         Smat[[k,q]]=Total[prob\[CurlyPhi]];
        ,Hmat[[k,q]]=Total[prob\[CurlyPhi]*(\[CapitalLambda]-E0)^(k+q-2)*\[CapitalLambda]];
         Smat[[k,q]]=Total[prob\[CurlyPhi]*(\[CapitalLambda]-E0)^(k+q-2)];];
   ,{q,1,d}];
,{k,1,d}];
Hmat=(Hmat+ConjugateTranspose[Hmat])/2.;
Smat=(Smat+ConjugateTranspose[Smat])/2.;
{Hmat,Smat}]


(*2. Chebyshev polynomial*)


funMatCP[\[CapitalLambda]_,E0_,d_,prob\[CurlyPhi]_,htot_]:=Module[{Hmat,Smat},
Hmat=ConstantArray[0,{d,d}];
Smat=ConstantArray[0,{d,d}];
Do[
   Do[
      Hmat[[k,q]]=Total[prob\[CurlyPhi]*ChebyshevT[k-1,(\[CapitalLambda]-E0)/htot]*ChebyshevT[q-1,(\[CapitalLambda]-E0)/htot]*\[CapitalLambda]];
      Smat[[k,q]]=Total[prob\[CurlyPhi]*ChebyshevT[k-1,(\[CapitalLambda]-E0)/htot]*ChebyshevT[q-1,(\[CapitalLambda]-E0)/htot]];
   ,{q,1,d}];
,{k,1,d}];
Hmat=(Hmat+ConjugateTranspose[Hmat])/2.;
Smat=(Smat+ConjugateTranspose[Smat])/2.;
{Hmat,Smat}]


(*3. Gaussian power*)


funMatGP[\[CapitalLambda]_,E0_,\[Tau]_,d_,prob\[CurlyPhi]_,htot_]:=Module[{Hmat,Smat,probG\[CurlyPhi],costList},
	costList=funCost[htot,\[Tau],d];
	probG\[CurlyPhi]=Exp[-(\[CapitalLambda]-E0)^2*\[Tau]^2]*prob\[CurlyPhi];
    Hmat=ConstantArray[0,{d,d}];
    Smat=ConstantArray[0,{d,d}];
	Do[
		Do[
			If[k+q-2==0,
			   Hmat[[k,q]]=Total[probG\[CurlyPhi]*\[CapitalLambda]]/(costList[[k]]*costList[[q]]);
			   Smat[[k,q]]=Total[probG\[CurlyPhi]]/(costList[[k]]*costList[[q]]);
			  ,Hmat[[k,q]]=Total[probG\[CurlyPhi]*(\[CapitalLambda]-E0)^(k+q-2)*\[CapitalLambda]]/(costList[[k]]*costList[[q]]);
			   Smat[[k,q]]=Total[probG\[CurlyPhi]*(\[CapitalLambda]-E0)^(k+q-2)]/(costList[[k]]*costList[[q]]);];
		,{q,1,d}]
	,{k,1,d}];
	Hmat=(Hmat+ConjugateTranspose[Hmat])/2.;
	Smat=(Smat+ConjugateTranspose[Smat])/2.;
{Hmat,Smat}]


(*funMatGP[\[CapitalLambda]_,E0_,\[Tau]_,d_,prob\[CurlyPhi]_]:=Module[{Hmat,Smat,probG\[CurlyPhi]},
probG\[CurlyPhi]=Exp[-(\[CapitalLambda]-E0)^2*\[Tau]^2]*prob\[CurlyPhi];
Hmat=ConstantArray[0,{d,d}];
Smat=ConstantArray[0,{d,d}];
Do[
   Do[
      If[k+q-2==0,
         Hmat[[k,q]]=Total[Re[probG\[CurlyPhi]*\[CapitalLambda]]]/(funCostGP[\[Tau],k-1]*funCostGP[\[Tau],q-1]);
         Smat[[k,q]]=Total[Re[probG\[CurlyPhi]]]/(funCostGP[\[Tau],k-1]*funCostGP[\[Tau],q-1]);
        ,Hmat[[k,q]]=Total[Re[probG\[CurlyPhi]*(\[CapitalLambda]-E0)^(k+q-2)*\[CapitalLambda]]]/(funCostGP[\[Tau],k-1]*funCostGP[\[Tau],q-1]);
         Smat[[k,q]]=Total[Re[probG\[CurlyPhi]*(\[CapitalLambda]-E0)^(k+q-2)]]/(funCostGP[\[Tau],k-1]*funCostGP[\[Tau],q-1]);
      ];
   ,{q,1,d}];
,{k,1,d}];
Hmat=(Hmat+ConjugateTranspose[Hmat])/2.;
Smat=(Smat+ConjugateTranspose[Smat])/2.;
{Hmat,Smat}
]*)


(*4. Inverse power*)


funMatIP[\[CapitalLambda]_,E0_,d_,prob\[CurlyPhi]_]:=Module[{Hmat,Smat},
Hmat=ConstantArray[0,{d,d}];
Smat=ConstantArray[0,{d,d}];
Do[
   Do[
      If[k+q-2==0,
         Hmat[[k,q]]=Total[prob\[CurlyPhi]*\[CapitalLambda]];
         Smat[[k,q]]=Total[prob\[CurlyPhi]];
        ,Hmat[[k,q]]=Total[prob\[CurlyPhi]*(\[CapitalLambda]-E0)^(-k-q+2)*\[CapitalLambda]];
         Smat[[k,q]]=Total[prob\[CurlyPhi]*(\[CapitalLambda]-E0)^(-k-q+2)];];
   ,{q,1,d}];
,{k,1,d}];
Hmat=(Hmat+ConjugateTranspose[Hmat])/2.;
Smat=(Smat+ConjugateTranspose[Smat])/2.;
{Hmat,Smat}
]


(*5. Imaginary time evolution*)


funMatITE[\[CapitalLambda]_,E0_,\[Tau]_,d_,prob\[CurlyPhi]_]:=Module[{Hmat,Smat,ITE},
ITE=Exp[-(\[CapitalLambda]-E0)*\[Tau]];
Hmat=ConstantArray[0,{d,d}];
Smat=ConstantArray[0,{d,d}];
Do[
   Do[
      Hmat[[k,q]]=Total[prob\[CurlyPhi]*ITE^(k+q-2)*\[CapitalLambda]];
      Smat[[k,q]]=Total[prob\[CurlyPhi]*ITE^(k+q-2)];
   ,{q,1,d}];
,{k,1,d}];
Hmat=(Hmat+ConjugateTranspose[Hmat])/2.;
Smat=(Smat+ConjugateTranspose[Smat])/2.;
{Hmat,Smat}]


(*6. Real time evolution*)


funMatRTE[\[CapitalLambda]_,E0_,\[CapitalDelta]t_,d_,prob\[CurlyPhi]_]:=Module[{Hmat,Smat,RTE},
RTE=Exp[I*(\[CapitalLambda]-E0)*\[CapitalDelta]t];
Hmat=ConstantArray[0,{d,d}];
Smat=ConstantArray[0,{d,d}];
Do[
   Do[
      Hmat[[k,q]]=Total[prob\[CurlyPhi]*RTE^(k-q)*\[CapitalLambda]];(*Hmat and Smat are complex Hermitian-Toeplitz matrices!!!*)
      Smat[[k,q]]=Total[prob\[CurlyPhi]*RTE^(k-q)];
   ,{q,1,d}];
,{k,1,d}];
Hmat=(Hmat+ConjugateTranspose[Hmat])/2.;
Smat=(Smat+ConjugateTranspose[Smat])/2.;
{Hmat,Smat}]


(*7. Filter*)


funMatF[\[CapitalLambda]_,E0_,\[CapitalDelta]E_,\[Tau]_,d_,prob\[CurlyPhi]_]:=Module[{Hmat,Smat},
Hmat=ConstantArray[0,{d,d}];
Smat=ConstantArray[0,{d,d}];
Do[
   Do[
      Hmat[[k,q]]=Total[prob\[CurlyPhi]*Sinc[(\[CapitalLambda]-(E0+(k-1)*\[CapitalDelta]E))*\[Tau]]*Sinc[(\[CapitalLambda]-(E0+(q-1)*\[CapitalDelta]E))*\[Tau]]*\[CapitalLambda]];
      Smat[[k,q]]=Total[prob\[CurlyPhi]*Sinc[(\[CapitalLambda]-(E0+(k-1)*\[CapitalDelta]E))*\[Tau]]*Sinc[(\[CapitalLambda]-(E0+(q-1)*\[CapitalDelta]E))*\[Tau]]];
   ,{q,1,d}];
,{k,1,d}];
Hmat=(Hmat+ConjugateTranspose[Hmat])/2.;
Smat=(Smat+ConjugateTranspose[Smat])/2.;
{Hmat,Smat}]


(*funTransformW[L_,d_,\[CapitalDelta]t_,E0_,\[CapitalDelta]E_]:=Module[{W},
W=ConstantArray[0,{L,d}];
Do[
   Do[
      W[[i,j]]=Exp[-I*(i-(L+1)/2)*\[CapitalDelta]t*(E0+\[CapitalDelta]E*(j-1))]/L;
   ,{j,1,d}];
,{i,1,L}];
W(*Dimension is L*d*)
]*)


(*funMatTAF[\[CapitalLambda]_,\[CapitalDelta]t_,L_,prob\[CurlyPhi]_,d_,E0_,\[CapitalDelta]E_]:=Module[{Hmat,Smat,W},
{Hmat,Smat}=funMatRTE[\[CapitalLambda],\[CapitalDelta]t,L,prob\[CurlyPhi]];
W=funTransformW[L,d,\[CapitalDelta]t,E0,\[CapitalDelta]E];
Hmat=ConjugateTranspose[W].Hmat.W;
Smat=ConjugateTranspose[W].Smat.W;
Hmat=(Hmat+ConjugateTranspose[Hmat])/2.;
Smat=(Smat+ConjugateTranspose[Smat])/2.;
{Hmat,Smat} (*Dimension is d*d, instead of L*L*)
]*)


(*Plot*)


funEpsilonGamma[Hmat_,Smat_,costH_,costS_,Id_,\[Eta]List_,Eg_,pg_]:=Module[{\[Epsilon]List,\[Gamma]List,\[Eta],EK,cn,\[Epsilon],\[Gamma]},
\[Epsilon]List=ConstantArray[0,Length[\[Eta]List]];
\[Gamma]List=ConstantArray[0,Length[\[Eta]List]];
Do[
\[Eta]=\[Eta]List[[j]];
{EK,cn}=funSubDiag[Hmat+2.*costH*\[Eta]*Id,Smat+2.*costS*\[Eta]*Id];
\[Epsilon]=EK-Eg;
\[Epsilon]List[[j]]=\[Epsilon];
\[Gamma]=(pg^2*\[Epsilon]^2)/(16*\[Eta]^2);
\[Gamma]List[[j]]=\[Gamma];
(*Print[{j,\[Eta],\[Gamma],\[Epsilon]}];*)
,{j,1,Length[\[Eta]List]}];
{\[Epsilon]List,\[Gamma]List}]


funEpsilonM[Hmat_,Smat_,costH_,costS_,Id_,\[Eta]List_,Eg_,d_,\[Kappa]_]:=Module[{\[Epsilon]List,MList,\[Eta],EK,cn,\[Epsilon]},
\[Epsilon]List=ConstantArray[0,Length[\[Eta]List]];
MList=ConstantArray[0,Length[\[Eta]List]];
Do[
\[Eta]=\[Eta]List[[j]];
MList[[j]]=0.5*d*Log[4d/\[Kappa]]/\[Eta]^2;
{EK,cn}=funSubDiag[Hmat+2.*costH*\[Eta]*Id,Smat+2.*costS*\[Eta]*Id];
\[Epsilon]=EK-Eg;
\[Epsilon]List[[j]]=\[Epsilon];
,{j,1,Length[\[Eta]List]}];
{MList,\[Epsilon]List}]


funEpsilonMRTE[Hmat_,Smat_,costH_,costS_,Id_,\[Eta]List_,Eg_,d_,\[Kappa]_]:=Module[{\[Epsilon]List,MList,\[Eta],EK,cn,\[Epsilon]},
\[Epsilon]List=ConstantArray[0,Length[\[Eta]List]];
MList=ConstantArray[0,Length[\[Eta]List]];
Do[
\[Eta]=\[Eta]List[[j]];
MList[[j]]=0.5*(2d-1)*Log[4d/\[Kappa]]/\[Eta]^2;
{EK,cn}=funSubDiag[Hmat+2.*costH*\[Eta]*Id,Smat+2.*costS*\[Eta]*Id];
\[Epsilon]=EK-Eg;
\[Epsilon]List[[j]]=\[Epsilon];
,{j,1,Length[\[Eta]List]}];
{MList,\[Epsilon]List}]


(*Regularisation in practice*)


funRegPrac[MList_,rep_,Hmat_,Smat_,d_,costH_,costS_,Id_,Eg_,pg_,complex_]:=Module[{\[Epsilon]List,\[Gamma]List,\[Sigma],hatH,hatS,\[Eta],hatEK,cn,\[Epsilon]},
\[Epsilon]List=ConstantArray[0,{rep,Length[MList]}];
\[Gamma]List=ConstantArray[0,{rep,Length[MList]}];
Do[
   Do[
      \[Sigma]=1/Sqrt[MList[[j]]];
      hatH=Hmat+RandomVariate[NormalDistribution[0,costH*\[Sigma]],{d,d}]+I*complex*RandomVariate[NormalDistribution[0,costH*\[Sigma]],{d,d}];
      hatH=(hatH+ConjugateTranspose[hatH])/2;(*Hermitian Gaussian noise matrix with no Hankel structure to stay up with references*)
      hatS=Smat+RandomVariate[NormalDistribution[0,costS*\[Sigma]],{d,d}]+I*complex*RandomVariate[NormalDistribution[0,costS*\[Sigma]],{d,d}];
      hatS=(hatS+ConjugateTranspose[hatS])/2;
      \[Eta]=Max[Norm[hatH-Hmat,2]/costH,Norm[hatS-Smat,2]/costS];(*optimal \[Eta]*)
      {hatEK,cn}=funSubDiag[hatH+costH*\[Eta]*Id,hatS+costS*\[Eta]*Id];
      \[Epsilon]=hatEK-Eg;
      \[Gamma]List[[i,j]]=(pg^2*\[Epsilon]^2)/(16*\[Eta]^2);
      \[Epsilon]List[[i,j]]=\[Epsilon];
   ,{j,1,Length[MList]}];
,{i,1,rep}];
{\[Epsilon]List,\[Gamma]List}]


(*Eta vs RMSE*)


funEtaGPRMSE[MList_,rep_,\[Eta]List_,Hmat_,Smat_,d_,costH_,costS_,Id_,Eg_]:=Module[{\[Epsilon]List,RMSEList,\[Sigma],hatH,hatS,\[Eta],hatEK,cn,\[Epsilon],GaussNoiList,HankelNoiMat},
RMSEList=ConstantArray[0,{Length[MList],Length[\[Eta]List]}];
Do[
   Do[
      \[Eta]=\[Eta]List[[j]];
      \[Epsilon]List=ConstantArray[0,rep];
      Do[
      \[Sigma]=1/Sqrt[MList[[i]]];      
      (*construct Hankel Gaussian Noise matrix of H*)
      GaussNoiList=RandomVariate[NormalDistribution[0,costH*\[Sigma]],2*d-1];
      HankelNoiMat=ConstantArray[0,{d,d}];
      Do[
         Do[
            HankelNoiMat[[i,j]]=GaussNoiList[[i+j-1]];
            ,{j,1,d}];
      ,{i,1,d}];
      hatH=Hmat+HankelNoiMat;
      (*construct Hankel Gaussian Noise matrix of S*)
      GaussNoiList=RandomVariate[NormalDistribution[0,costS*\[Sigma]],2*d-1];
      HankelNoiMat=ConstantArray[0,{d,d}];
      Do[
         Do[
            HankelNoiMat[[i,j]]=GaussNoiList[[i+j-1]];
            ,{j,1,d}];
      ,{i,1,d}];
      hatS=Smat+HankelNoiMat;
      {hatEK,cn}=funSubDiag[hatH+costH*\[Eta]*Id,hatS+costS*\[Eta]*Id];
      \[Epsilon]=hatEK-Eg;
      \[Epsilon]List[[k]]=\[Epsilon];
      ,{k,1,rep}];
      RMSEList[[i,j]]=RootMeanSquare[\[Epsilon]List];
   ,{j,1,Length[\[Eta]List]}];
,{i,1,Length[MList]}];
RMSEList]


funEtaRTERMSE[MList_,rep_,\[Eta]List_,Hmat_,Smat_,d_,costH_,costS_,Id_,Eg_]:=Module[{\[Epsilon]List,RMSEList,\[Sigma],hatH,hatS,\[Eta],hatEK,cn,\[Epsilon],GaussNoiList,ToeplitzNoiMat},
RMSEList=ConstantArray[0,{Length[MList],Length[\[Eta]List]}];
Do[
   Do[
      \[Eta]=\[Eta]List[[j]];
      \[Epsilon]List=ConstantArray[0,rep];
      Do[
      \[Sigma]=1/Sqrt[MList[[i]]];      
      (*construct complex Hermite-Toeplitz Gaussian Noise matrix of H*)
      GaussNoiList=RandomVariate[NormalDistribution[0,costH*\[Sigma]],2*d-1];
      ToeplitzNoiMat=ConstantArray[0,{d,d}];
       Do[
         Do[
            ToeplitzNoiMat[[i,j]]=GaussNoiList[[Abs[i-j]+1]]+Sign[i-j]*I*GaussNoiList[[Abs[i-j]+d]];
            ,{j,1,d}];
      ,{i,1,d}];
      hatH=Hmat+ToeplitzNoiMat;
      (*construct complex Hermite-Toeplitz Gaussian Noise matrix of S*)
      GaussNoiList=RandomVariate[NormalDistribution[0,costS*\[Sigma]],2*d-1];
      ToeplitzNoiMat=ConstantArray[0,{d,d}];
       Do[
         Do[
            ToeplitzNoiMat[[i,j]]=GaussNoiList[[Abs[i-j]+1]]+Sign[i-j]*I*GaussNoiList[[Abs[i-j]+d]];
            ,{j,1,d}];
      ,{i,1,d}];
      hatS=Smat+ToeplitzNoiMat;
      {hatEK,cn}=funSubDiag[hatH+costH*\[Eta]*Id,hatS+costS*\[Eta]*Id];
      \[Epsilon]=hatEK-Eg;
      \[Epsilon]List[[k]]=\[Epsilon];
      ,{k,1,rep}];
      RMSEList[[i,j]]=RootMeanSquare[\[Epsilon]List];
   ,{j,1,Length[\[Eta]List]}];
,{i,1,Length[MList]}];
RMSEList]


funEtaFRMSE[MList_,rep_,\[Eta]List_,Hmat_,Smat_,d_,costH_,costS_,Id_,Eg_]:=Module[{\[Epsilon]List,RMSEList,\[Sigma],hatH,hatS,\[Eta],hatEK,cn,\[Epsilon]},
RMSEList=ConstantArray[0,{Length[MList],Length[\[Eta]List]}];
Do[
   Do[
      \[Eta]=\[Eta]List[[j]];
      \[Epsilon]List=ConstantArray[0,rep];
      Do[
      \[Sigma]=1/Sqrt[MList[[i]]];      
      hatH=Hmat+RandomVariate[NormalDistribution[0,costH*\[Sigma]],{d,d}];
      hatH=(hatH+ConjugateTranspose[hatH])/2;(*Real Hermitian Gaussian noise matrix*)
      hatS=Smat+RandomVariate[NormalDistribution[0,costS*\[Sigma]],{d,d}];
      hatS=(hatS+ConjugateTranspose[hatS])/2;
      {hatEK,cn}=funSubDiag[hatH+costH*\[Eta]*Id,hatS+costS*\[Eta]*Id];
      \[Epsilon]=hatEK-Eg;
      \[Epsilon]List[[k]]=\[Epsilon];
      ,{k,1,rep}];
      RMSEList[[i,j]]=RootMeanSquare[\[Epsilon]List];
   ,{j,1,Length[\[Eta]List]}];
,{i,1,Length[MList]}];
RMSEList]


(*eta vs Abs[epsilon] in 1-kappa probability*)


(**suitable for P, GP, IP and ITE*)
funEtaEpsilonGP[MList_,rep_,\[Kappa]_,\[Eta]List_,Hmat_,Smat_,d_,costH_,costS_,Id_,Eg_]:=Module[{\[Epsilon]List,EpsilonList,\[Sigma],M,hatH,hatS,\[Eta],\[Eta]pracList,hatEK,cn,\[Epsilon],GaussNoiList,HankelNoiMat},
EpsilonList=ConstantArray[0,{Length[MList],Length[\[Eta]List]}];
\[Eta]pracList=ConstantArray[0,Length[MList]];
Do[
M=MList[[i]];
\[Sigma]=1/Sqrt[M]; 
hatH=ConstantArray[0,{rep,d,d}];
hatS=ConstantArray[0,{rep,d,d}];   
   Do[
      (*construct Hankel Gaussian Noise matrix of H*)
      GaussNoiList=RandomVariate[NormalDistribution[0,costH*\[Sigma]],2*d-1];
      HankelNoiMat=ConstantArray[0,{d,d}];
      Do[
         Do[
            HankelNoiMat[[i,j]]=GaussNoiList[[i+j-1]];
            ,{j,1,d}];
      ,{i,1,d}];
      hatH[[k]]=Hmat+HankelNoiMat;
      (*construct Hankel Gaussian Noise matrix of S*)
      GaussNoiList=RandomVariate[NormalDistribution[0,costS*\[Sigma]],2*d-1];
      HankelNoiMat=ConstantArray[0,{d,d}];
      Do[
         Do[
            HankelNoiMat[[i,j]]=GaussNoiList[[i+j-1]];
            ,{j,1,d}];
      ,{i,1,d}];
      hatS[[k]]=Smat+HankelNoiMat;
   ,{k,1,rep}];
   Do[
      \[Eta]=\[Eta]List[[j]];
      \[Epsilon]List=ConstantArray[0,rep];
      Do[     
      {hatEK,cn}=funSubDiag[hatH[[k]]+costH*\[Eta]*Id,hatS[[k]]+costS*\[Eta]*Id];
      \[Epsilon]=Abs[hatEK-Eg];
      \[Epsilon]List[[k]]=\[Epsilon];
      ,{k,1,rep}];
      EpsilonList[[i,j]]=Sort[\[Epsilon]List,Less][[Round[(1-\[Kappa])*rep]]];
   ,{j,1,Length[\[Eta]List]}];
\[Eta]pracList[[i]]=Sqrt[0.5*d/M Log[4*d/\[Kappa]]];(*Matrix Gaussian series for real Hermite/Hankel matrix*)
,{i,1,Length[MList]}];
{EpsilonList,\[Eta]pracList}]


funEtaEpsilonRTE[MList_,rep_,\[Kappa]_,\[Eta]List_,Hmat_,Smat_,d_,costH_,costS_,Id_,Eg_]:=Module[{\[Epsilon]List,EpsilonList,\[Sigma],M,hatH,hatS,\[Eta],\[Eta]pracList,hatEK,cn,\[Epsilon],GaussNoiList,ToeplitzNoiMat},
EpsilonList=ConstantArray[0,{Length[MList],Length[\[Eta]List]}];
\[Eta]pracList=ConstantArray[0,Length[MList]];
Do[
M=MList[[i]];
\[Sigma]=1/Sqrt[M];  
hatH=ConstantArray[0,{rep,d,d}];
hatS=ConstantArray[0,{rep,d,d}];   
   Do[    
      (*construct complex Hermite-Toeplitz Gaussian Noise matrix of H*)
      GaussNoiList=RandomVariate[NormalDistribution[0,costH*\[Sigma]],2*d-1];
      ToeplitzNoiMat=ConstantArray[0,{d,d}];
       Do[
         Do[
            ToeplitzNoiMat[[i,j]]=GaussNoiList[[Abs[i-j]+1]]+Sign[i-j]*I*GaussNoiList[[Abs[i-j]+d]];
            ,{j,1,d}];
      ,{i,1,d}];
      hatH[[k]]=Hmat+ToeplitzNoiMat;
      (*construct complex Hermite-Toeplitz Gaussian Noise matrix of S*)
      GaussNoiList=RandomVariate[NormalDistribution[0,costS*\[Sigma]],2*d-1];
      ToeplitzNoiMat=ConstantArray[0,{d,d}];
       Do[
         Do[
            ToeplitzNoiMat[[i,j]]=GaussNoiList[[Abs[i-j]+1]]+Sign[i-j]*I*GaussNoiList[[Abs[i-j]+d]];
            ,{j,1,d}];
      ,{i,1,d}];
      hatS[[k]]=Smat+ToeplitzNoiMat;
   ,{k,1,rep}];
   Do[
      \[Eta]=\[Eta]List[[j]];
      \[Epsilon]List=ConstantArray[0,rep];
      Do[
      {hatEK,cn}=funSubDiag[hatH[[k]]+costH*\[Eta]*Id,hatS[[k]]+costS*\[Eta]*Id];
      \[Epsilon]=Abs[hatEK-Eg];
      \[Epsilon]List[[k]]=\[Epsilon];
      ,{k,1,rep}];
      EpsilonList[[i,j]]=Sort[\[Epsilon]List,Less][[Round[(1-\[Kappa])*rep]]];
   ,{j,1,Length[\[Eta]List]}];
\[Eta]pracList[[i]]=Sqrt[0.5*(2d-1)/M Log[4*d/\[Kappa]]];(*Matrix Gaussian series for complex Hermite/Toeplitz matrix*)
,{i,1,Length[MList]}];
{EpsilonList,\[Eta]pracList}]


(*suitable for CP and F*)
funEtaEpsilonF[MList_,rep_,\[Kappa]_,\[Eta]List_,Hmat_,Smat_,d_,costH_,costS_,Id_,Eg_]:=Module[{\[Epsilon]List,EpsilonList,\[Sigma],M,hatH,hatS,\[Eta],\[Eta]pracList,hatEK,cn,\[Epsilon]},
EpsilonList=ConstantArray[0,{Length[MList],Length[\[Eta]List]}];
\[Eta]pracList=ConstantArray[0,Length[MList]];
Do[
M=MList[[i]];
\[Sigma]=1/Sqrt[M];  
hatH=ConstantArray[0,{rep,d,d}];
hatS=ConstantArray[0,{rep,d,d}];   
   Do[
      hatH[[k]]=Hmat+RandomVariate[NormalDistribution[0,costH*\[Sigma]],{d,d}];
      hatH[[k]]=(hatH[[k]]+ConjugateTranspose[hatH[[k]]])/2;(*Real Hermitian Gaussian noise matrix*)
      hatS[[k]]=Smat+RandomVariate[NormalDistribution[0,costS*\[Sigma]],{d,d}];
      hatS[[k]]=(hatS[[k]]+ConjugateTranspose[hatS[[k]]])/2;
   ,{k,1,rep}];
   Do[
      \[Eta]=\[Eta]List[[j]];
      \[Epsilon]List=ConstantArray[0,rep];
      Do[    
      {hatEK,cn}=funSubDiag[hatH[[k]]+costH*\[Eta]*Id,hatS[[k]]+costS*\[Eta]*Id];
      \[Epsilon]=hatEK-Eg;
      \[Epsilon]List[[k]]=\[Epsilon];
      ,{k,1,rep}];
      EpsilonList[[i,j]]=Sort[\[Epsilon]List,Less][[Round[(1-\[Kappa])*rep]]];
   ,{j,1,Length[\[Eta]List]}];
\[Eta]pracList[[i]]=Sqrt[0.5*d/M Log[4*d/\[Kappa]]];(*Matrix Gaussian series for real Hermite/Hankel matrix*)
,{i,1,Length[MList]}];
{EpsilonList,\[Eta]pracList}]


(*Practical eta*)(*M vs epsilon in 1-kappa probability using practical eta*)


funExtract[\[Epsilon]List_,MList_,rep_,\[Kappa]_]:=Module[{\[Epsilon]List\[Kappa]},
\[Epsilon]List\[Kappa]=ConstantArray[0,Length[MList]];
Do[
\[Epsilon]List\[Kappa][[i]]=Sort[\[Epsilon]List[[;;,i]],Less][[Round[(1-\[Kappa])*rep]]];
,{i,1,Length[MList]}];
\[Epsilon]List\[Kappa]]


(*GP,P,IP,ITE*)
funEtaPracGP[MList_,rep_,n_,Hmat_,Smat_,d_,\[Kappa]_,costH_,costS_,Id_,Eg_,pg_]:=Module[{\[Epsilon]List,\[Gamma]List,M,\[Sigma],hatH,hatS,\[Eta]prac,\[Eta],hatEK,cn,\[Epsilon],GaussNoiList,HankelNoiMat},
\[Epsilon]List=ConstantArray[0,{rep,Length[MList]}];
\[Gamma]List=ConstantArray[0,{rep,Length[MList]}];
Do[
   Do[
      M=MList[[j]];
      \[Sigma]=1/Sqrt[M];
      (*construct Hankel Gaussian Noise matrix of H*)
      GaussNoiList=RandomVariate[NormalDistribution[0,costH*\[Sigma]],2*d-1];
      HankelNoiMat=ConstantArray[0,{d,d}];
      Do[
         Do[
            HankelNoiMat[[i,j]]=GaussNoiList[[i+j-1]];
            ,{j,1,d}];
      ,{i,1,d}];
      hatH=Hmat+HankelNoiMat;
      (*construct Hankel Gaussian Noise matrix of S*)
      GaussNoiList=RandomVariate[NormalDistribution[0,costS*\[Sigma]],2*d-1];
      HankelNoiMat=ConstantArray[0,{d,d}];
      Do[
         Do[
            HankelNoiMat[[i,j]]=GaussNoiList[[i+j-1]];
            ,{j,1,d}];
      ,{i,1,d}];
      hatS=Smat+HankelNoiMat;
      \[Eta]prac=Sqrt[0.5*d/M Log[4*d/\[Kappa]]];(*Matrix Gaussian series for real Hermite/Hankel matrix*)
      \[Eta]=n*\[Eta]prac;
      {hatEK,cn}=funSubDiag[hatH+costH*\[Eta]*Id,hatS+costS*\[Eta]*Id];
      \[Epsilon]=hatEK-Eg;
      \[Gamma]List[[i,j]]=(pg^2*\[Epsilon]^2)/(16*\[Eta]^2);
      \[Epsilon]List[[i,j]]=\[Epsilon];
   ,{j,1,Length[MList]}];
,{i,1,rep}];
{\[Epsilon]List,\[Gamma]List}]


funEtaPracRTE[MList_,rep_,n_,Hmat_,Smat_,d_,\[Kappa]_,costH_,costS_,Id_,Eg_,pg_]:=Module[{\[Epsilon]List,\[Gamma]List,M,\[Sigma],hatH,hatS,\[Eta]prac,\[Eta],hatEK,cn,\[Epsilon],GaussNoiList,ToeplitzNoiMat},
\[Epsilon]List=ConstantArray[0,{rep,Length[MList]}];
\[Gamma]List=ConstantArray[0,{rep,Length[MList]}];
Do[
   Do[
      M=MList[[j]];
      \[Sigma]=1/Sqrt[M];
      (*construct complex Hermite-Toeplitz Gaussian Noise matrix of H*)
      GaussNoiList=RandomVariate[NormalDistribution[0,costH*\[Sigma]],2*d-1];
      ToeplitzNoiMat=ConstantArray[0,{d,d}];
       Do[
         Do[
            ToeplitzNoiMat[[i,j]]=GaussNoiList[[Abs[i-j]+1]]+Sign[i-j]*I*GaussNoiList[[Abs[i-j]+d]];
            ,{j,1,d}];
      ,{i,1,d}];
      hatH=Hmat+ToeplitzNoiMat;
      (*construct complex Hermite-Toeplitz Gaussian Noise matrix of S*)
      GaussNoiList=RandomVariate[NormalDistribution[0,costS*\[Sigma]],2*d-1];
      ToeplitzNoiMat=ConstantArray[0,{d,d}];
       Do[
         Do[
            ToeplitzNoiMat[[i,j]]=GaussNoiList[[Abs[i-j]+1]]+Sign[i-j]*I*GaussNoiList[[Abs[i-j]+d]];
            ,{j,1,d}];
      ,{i,1,d}];
      hatS=Smat+ToeplitzNoiMat;
      \[Eta]prac=Sqrt[0.5*(2d-1)/M Log[4*d/\[Kappa]]];(*Matrix Gaussian series for complex Hermite/Toeplitz matrix*)
      \[Eta]=n*\[Eta]prac;
      {hatEK,cn}=funSubDiag[hatH+costH*\[Eta]*Id,hatS+costS*\[Eta]*Id];
      \[Epsilon]=hatEK-Eg;
      \[Gamma]List[[i,j]]=(pg^2*\[Epsilon]^2)/(16*\[Eta]^2);
      \[Epsilon]List[[i,j]]=\[Epsilon];
   ,{j,1,Length[MList]}];
,{i,1,rep}];
{\[Epsilon]List,\[Gamma]List}]


(*F,CP,rescaled GP*)
funEtaPracF[MList_,rep_,n_,Hmat_,Smat_,d_,\[Kappa]_,costH_,costS_,Id_,Eg_,pg_]:=Module[{\[Epsilon]List,\[Gamma]List,M,\[Sigma],hatH,hatS,\[Eta]prac,\[Eta],hatEK,cn,\[Epsilon]},
\[Epsilon]List=ConstantArray[0,{rep,Length[MList]}];
\[Gamma]List=ConstantArray[0,{rep,Length[MList]}];
Do[
   Do[
      M=MList[[j]];
      \[Sigma]=1/Sqrt[M];
      hatH=Hmat+RandomVariate[NormalDistribution[0,costH*\[Sigma]],{d,d}];
      hatH=(hatH+ConjugateTranspose[hatH])/2;(*Real Hermitian Gaussian noise matrix*)
      hatS=Smat+RandomVariate[NormalDistribution[0,costS*\[Sigma]],{d,d}];
      hatS=(hatS+ConjugateTranspose[hatS])/2;
      \[Eta]prac=Sqrt[0.5*d/M Log[4*d/\[Kappa]]];(*Matrix Gaussian series for real Hermite/Hankel matrix*)
      \[Eta]=n*\[Eta]prac;
      {hatEK,cn}=funSubDiag[hatH+costH*\[Eta]*Id,hatS+costS*\[Eta]*Id];
      \[Epsilon]=hatEK-Eg;
      \[Gamma]List[[i,j]]=(pg^2*\[Epsilon]^2)/(16*\[Eta]^2);
      \[Epsilon]List[[i,j]]=\[Epsilon];
   ,{j,1,Length[MList]}];
,{i,1,rep}];
{\[Epsilon]List,\[Gamma]List}]


(*Thresholding*)


thresholding[hatH_,hatS_,\[Theta]_]:=Module[{vals,vecs,index,V,Hth,Sth,hatEK,cn},
{vals,vecs}=funSpectrum[hatS];
index={};
Do[
If[vals[[k]]>\[Theta],AppendTo[index,k]];
,{k,1,d}];
V=Transpose[Extract[vecs,Transpose[{index}]]];
Hth=ConjugateTranspose[V] . hatH . V;
Hth=(Hth+ConjugateTranspose[Hth])/2;
Sth=ConjugateTranspose[V] . hatS . V;
Sth=(Sth+ConjugateTranspose[Sth])/2;
{hatEK,cn}=funSubDiag[Hth,Sth];
hatEK]


funThrPracGP[MList_,rep_,Hmat_,Smat_,d_,costH_,costS_,th_,Eg_]:=Module[{\[Epsilon]List,\[Sigma],hatH,hatS,\[Theta],hatEK,\[Epsilon],GaussNoiList,HankelNoiMat},
\[Epsilon]List=ConstantArray[0,{rep,Length[MList]}];
Do[
   Do[
      \[Sigma]=1/Sqrt[MList[[j]]];
      (*construct Hankel Gaussian Noise matrix of H*)
      GaussNoiList=RandomVariate[NormalDistribution[0,costH*\[Sigma]],2*d-1];
      HankelNoiMat=ConstantArray[0,{d,d}];
      Do[
         Do[
            HankelNoiMat[[i,j]]=GaussNoiList[[i+j-1]];
            ,{j,1,d}];
      ,{i,1,d}];
      hatH=Hmat+HankelNoiMat;
      (*construct Hankel Gaussian Noise matrix of S*)
      GaussNoiList=RandomVariate[NormalDistribution[0,costS*\[Sigma]],2*d-1];
      HankelNoiMat=ConstantArray[0,{d,d}];
      Do[
         Do[
            HankelNoiMat[[i,j]]=GaussNoiList[[i+j-1]];
            ,{j,1,d}];
      ,{i,1,d}];
      hatS=Smat+HankelNoiMat;
      \[Theta]=th*\[Sigma];
      hatEK=thresholding[hatH,hatS,\[Theta]];
      \[Epsilon]=Abs[hatEK-Eg];
      \[Epsilon]List[[i,j]]=\[Epsilon];
   ,{j,1,Length[MList]}];
,{i,1,rep}];
\[Epsilon]List]


funThrPracRTE[MList_,rep_,Hmat_,Smat_,d_,costH_,costS_,th_,Eg_]:=Module[{\[Epsilon]List,\[Sigma],hatH,hatS,\[Theta],hatEK,\[Epsilon],GaussNoiList,ToeplitzNoiMat},
\[Epsilon]List=ConstantArray[0,{rep,Length[MList]}];
Do[
   Do[
      \[Sigma]=1/Sqrt[MList[[j]]];
      (*construct complex Hermite-Toeplitz Gaussian Noise matrix of H*)
      GaussNoiList=RandomVariate[NormalDistribution[0,costH*\[Sigma]],2*d-1];
      ToeplitzNoiMat=ConstantArray[0,{d,d}];
       Do[
         Do[
            ToeplitzNoiMat[[i,j]]=GaussNoiList[[Abs[i-j]+1]]+Sign[i-j]*I*GaussNoiList[[Abs[i-j]+d]];
            ,{j,1,d}];
      ,{i,1,d}];
      hatH=Hmat+ToeplitzNoiMat;
      (*construct complex Hermite-Toeplitz Gaussian Noise matrix of S*)
      GaussNoiList=RandomVariate[NormalDistribution[0,costS*\[Sigma]],2*d-1];
      ToeplitzNoiMat=ConstantArray[0,{d,d}];
       Do[
         Do[
            ToeplitzNoiMat[[i,j]]=GaussNoiList[[Abs[i-j]+1]]+Sign[i-j]*I*GaussNoiList[[Abs[i-j]+d]];
            ,{j,1,d}];
      ,{i,1,d}];
      hatS=Smat+ToeplitzNoiMat;
      \[Theta]=th*\[Sigma];
      hatEK=thresholding[hatH,hatS,\[Theta]];
      \[Epsilon]=Abs[hatEK-Eg];
      \[Epsilon]List[[i,j]]=\[Epsilon];
   ,{j,1,Length[MList]}];
,{i,1,rep}];
\[Epsilon]List]


funThrPracF[MList_,rep_,Hmat_,Smat_,d_,costH_,costS_,th_,Eg_]:=Module[{\[Epsilon]List,\[Sigma],hatH,hatS,\[Theta],hatEK,\[Epsilon]},
\[Epsilon]List=ConstantArray[0,{rep,Length[MList]}];
Do[
   Do[
      \[Sigma]=1/Sqrt[MList[[j]]];
      hatH=Hmat+RandomVariate[NormalDistribution[0,costH*\[Sigma]],{d,d}];
      hatH=(hatH+ConjugateTranspose[hatH])/2;(*Real Hermitian Gaussian noise matrix*)
      hatS=Smat+RandomVariate[NormalDistribution[0,costS*\[Sigma]],{d,d}];
      hatS=(hatS+ConjugateTranspose[hatS])/2;
      \[Theta]=th*\[Sigma];
      hatEK=thresholding[hatH,hatS,\[Theta]];
      \[Epsilon]=Abs[hatEK-Eg];
      \[Epsilon]List[[i,j]]=\[Epsilon];
   ,{j,1,Length[MList]}];
,{i,1,rep}];
\[Epsilon]List]
