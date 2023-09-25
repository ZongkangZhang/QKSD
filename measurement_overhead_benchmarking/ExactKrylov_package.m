(* ::Package:: *)

(*Model*)


funGraph[Nv_,GraphType_]:=Module[{EL,vv},(
	If[StringMatchQ[GraphType,"Chain"],(
		EL={};
		Do[(
			AppendTo[EL,UndirectedEdge[v,Mod[v+1,Nv]]]
		),{v,0,Nv-2}];
	)];
	
	If[StringMatchQ[GraphType,"Ladder"],(
		EL={};
		Do[(
			AppendTo[EL,UndirectedEdge[v,Mod[v+1,Nv]]];
		),{v,0,Nv-1,2}];
		Do[(
			AppendTo[EL,UndirectedEdge[v,Mod[v+2,Nv]]];
			AppendTo[EL,UndirectedEdge[Mod[v+1,Nv],Mod[v+3,Nv]]]
		),{v,0,Nv-3,2}];
	)];
	
	If[!StringMatchQ[GraphType,"Chain"]&&!StringMatchQ[GraphType,"Ladder"],(
		EL={};
		Do[(
			Do[(
				vv=RandomChoice[Drop[Table[vv,{vv,0,Nv-1}],{v+1}]];
				AppendTo[EL,UndirectedEdge[v,vv]]
			),{i,1,ToExpression[GraphType]}];
		),{v,0,Nv-1}]
		(*RG=RandomGraph[{Nv,ToExpression[GraphType]*Nv}];
		EL=EdgeList[RG]*)
	)];
	
	Return[EL]
)]


funHeisenberg[Nq_,EL_]:=Module[{Istr,Model,q,qq,ps},(
	Istr="";
	Do[(
		Istr=StringInsert[Istr,"I",q+1];
	),{q,0,Nq-1}];
	Model={};
	
	Do[(
		q=Mod[EL[[l,1]],Nq];
		qq=Mod[EL[[l,2]],Nq];
		
		ps=StringReplacePart[Istr,"X",{q+1,q+1}];
		ps=StringReplacePart[ps,"X",{qq+1,qq+1}];
		AppendTo[Model,{1.,ps}];
		
		ps=StringReplacePart[Istr,"Y",{q+1,q+1}];
		ps=StringReplacePart[ps,"Y",{qq+1,qq+1}];
		AppendTo[Model,{1.,ps}];
		
		ps=StringReplacePart[Istr,"Z",{q+1,q+1}];
		ps=StringReplacePart[ps,"Z",{qq+1,qq+1}];
		AppendTo[Model,{1.,ps}];
	),{l,1,Length[EL]}];
	
	Return[Model]
)]


funFermiHubbard[Nq_,EL_,u_]:=Module[{Istr,Model,v,vv,q,qq,ps},(
	Istr="";
	Do[(
		Istr=StringInsert[Istr,"I",q+1];
	),{q,0,Nq-1}];
	Model={};
	
	Do[(
		v=Mod[EL[[l,1]],Nq/2];
		vv=Mod[EL[[l,2]],Nq/2];
		
		q=2*Min[{v,vv}];
		qq=2*Max[{v,vv}];
		ps=StringReplacePart[Istr,"X",{q+1,q+1}];
		Do[(
			ps=StringReplacePart[ps,"Z",{qqq+1,qqq+1}];
		),{qqq,q+1,qq-1}];
		ps=StringReplacePart[ps,"X",{qq+1,qq+1}];
		AppendTo[Model,{-1./2.,ps}];
		ps=StringReplacePart[Istr,"Y",{q+1,q+1}];
		Do[(
			ps=StringReplacePart[ps,"Z",{qqq+1,qqq+1}];
		),{qqq,q+1,qq-1}];
		ps=StringReplacePart[ps,"Y",{qq+1,qq+1}];
		AppendTo[Model,{-1./2.,ps}];
		
		q=2*Min[{v,vv}]+1;
		qq=2*Max[{v,vv}]+1;
		ps=StringReplacePart[Istr,"X",{q+1,q+1}];
		Do[(
			ps=StringReplacePart[ps,"Z",{qqq+1,qqq+1}];
		),{qqq,q+1,qq-1}];
		ps=StringReplacePart[ps,"X",{qq+1,qq+1}];
		AppendTo[Model,{-1./2.,ps}];
		ps=StringReplacePart[Istr,"Y",{q+1,q+1}];
		Do[(
			ps=StringReplacePart[ps,"Z",{qqq+1,qqq+1}];
		),{qqq,q+1,qq-1}];
		ps=StringReplacePart[ps,"Y",{qq+1,qq+1}];
		AppendTo[Model,{-1./2.,ps}];
		
		(*ps=StringReplacePart[Istr,"Z",{q+1,q+1}];
		ps=StringReplacePart[ps,"Z",{qq+1,qq+1}];
		AppendTo[Model,{1.,ps}];*)
	),{l,1,Length[EL]}];
	
	Do[(
		q=2*v;
		qq=2*v+1;
		ps=StringReplacePart[Istr,"Z",{q+1,q+1}];
		ps=StringReplacePart[ps,"Z",{qq+1,qq+1}];
		AppendTo[Model,{u/4.,ps}];
	),{v,0,Nq/2-1}];
	
	Return[Model]
)]


(*funFermiHubbard[Nq_,EL_,u_]:=Module[{Istr,Model,v,vv,q,qq,ps},(
	Istr="";
	Do[(
		Istr=StringInsert[Istr,"I",q+1];
	),{q,0,Nq-1}];
	Model={};
	
	Do[(
		v=Mod[EL[[l,1]],Nq/2];
		vv=Mod[EL[[l,2]],Nq/2];
		
		q=2*Min[{v,vv}];
		qq=2*Max[{v,vv}];
		ps=StringReplacePart[Istr,"X",{q+1,q+1}];
		Do[(
			ps=StringReplacePart[Istr,"Z",{q+1,q+1}];
		),{qqq,q+1,qq-1}];
		ps=StringReplacePart[ps,"X",{qq+1,qq+1}];
		AppendTo[Model,{-1./2.,ps}];
		ps=StringReplacePart[Istr,"Y",{q+1,q+1}];
		Do[(
			ps=StringReplacePart[Istr,"Z",{q+1,q+1}];
		),{qqq,q+1,qq-1}];
		ps=StringReplacePart[ps,"Y",{qq+1,qq+1}];
		AppendTo[Model,{-1./2.,ps}];
		
		q=2*Min[{v,vv}]+1;
		qq=2*Max[{v,vv}]+1;
		ps=StringReplacePart[Istr,"X",{q+1,q+1}];
		Do[(
			ps=StringReplacePart[Istr,"Z",{q+1,q+1}];
		),{qqq,q+1,qq-1}];
		ps=StringReplacePart[ps,"X",{qq+1,qq+1}];
		AppendTo[Model,{-1./2.,ps}];
		ps=StringReplacePart[Istr,"Y",{q+1,q+1}];
		Do[(
			ps=StringReplacePart[Istr,"Z",{q+1,q+1}];
		),{qqq,q+1,qq-1}];
		ps=StringReplacePart[ps,"Y",{qq+1,qq+1}];
		AppendTo[Model,{-1./2.,ps}];
		
		ps=StringReplacePart[Istr,"Z",{q+1,q+1}];
		ps=StringReplacePart[ps,"Z",{qq+1,qq+1}];
		AppendTo[Model,{1.,ps}];
	),{l,1,Length[EL]}];
	
	Do[(
		q=2*v;
		qq=2*v+1;
		ps=StringReplacePart[Istr,"Z",{q+1,q+1}];
		ps=StringReplacePart[ps,"Z",{qq+1,qq+1}];
		AppendTo[Model,{u/4.,ps}];
	),{v,0,Nq/2-1}];
	
	Return[Model]
)]*)


(*Reference state*)


funPairwiseSinglet[Nq_]:=Module[{U,qq,\[Psi]},(
	U=IdentityMatrix[2^Nq];
	Do[(
		qq=Mod[q+1,Nq];
		U=U . ((funPX[q]-funPX[qq])/N[Sqrt[2]]);
	),{q,0,Nq-1,2}];
	
	\[Psi]=Table[{0.},2^Nq];
	\[Psi][[1,1]]=1.;
	\[Psi]=U . \[Psi];
	
	Return[\[Psi]]
)]


funHartreeFock[Nq_,EL_]:=Module[{Model,Ham,EE,ES,\[Psi]},(
	Model=funFermiHubbard[Nq,EL,0.];
	Ham=funHamiltonianQubit[Model];
	{EE,ES}=funSpectrum[Ham];
	\[Psi]=Transpose[{ES[[1]]}];
	Return[\[Psi]]
)]


(*Diagonalisation*)


funDiagonalisation[Hmat_,Smat_]:=Module[{Svals,Svecs,cn,V,Heff,Hvals,Hvecs,EK,SK},(
	{Svals,Svecs}=funSpectrum[Smat];
	cn=Max[Abs[Svals]]/Min[Abs[Svals]];
	V=Transpose[Svecs] . DiagonalMatrix[1./Sqrt[Svals]];
	Heff=ConjugateTranspose[V] . Hmat . V;
	{Hvals,Hvecs}=funSpectrum[Heff];
	EK=Hvals[[1]];
	Return[{EK,cn}]
)]


(*Functions*)


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


(*funCost[\[Tau]_,k_]:=Module[{\[Chi],A,Cost},(
	\[Chi]=0.125;
	A=1.8946081370976193;
	If[k==1,(
		Cost=Sqrt[1./(1.-2.*\[Chi])];
	),(
		Cost=A*(k-1)^((k-1)/2)/Exp[(k-1)/2]/\[Tau]^(k-1);
	)];
	Return[Cost]
)]*)


(*funCost[\[Tau]_,k_]:=Module[{\[Chi],Cost},(
	\[Chi]=0.125;
	If[k==1,(
		Cost=Sqrt[1./(1.-2.*\[Chi])];
	),(
		Cost=Sqrt[2.*(k-1)!/(1.-4.*\[Chi])]/\[Tau]^(k-1);
	)];
	Return[Cost]
)]*)


(*Matrices*)


funMatPower[EE_,Pro\[Psi]_,d_,E0_]:=Module[{Hmat,Smat},(
	Hmat=Table[Table[0.,d],d];
	Smat=Table[Table[0.,d],d];
	Do[(
		Do[(
			If[k+q-2==0,(
				Hmat[[k,q]]=Total[Pro\[Psi]*EE];
				Smat[[k,q]]=Total[Pro\[Psi]]
			),(
				Hmat[[k,q]]=Total[Pro\[Psi]*(EE-E0)^(k+q-2)*EE];
				Smat[[k,q]]=Total[Pro\[Psi]*(EE-E0)^(k+q-2)]
			)]
			
		),{q,1,d}]
	),{k,1,d}];
	Hmat=(Hmat+ConjugateTranspose[Hmat])/2.;
	Smat=(Smat+ConjugateTranspose[Smat])/2.;
	Return[{Hmat,Smat}]
)]


funMatChebyshev[EE_,Pro\[Psi]_,d_,htot_,E0_]:=Module[{Hmat,Smat},(
	Hmat=Table[Table[0.,d],d];
	Smat=Table[Table[0.,d],d];
	Do[(
		Do[(
			Hmat[[k,q]]=Total[Pro\[Psi]*ChebyshevT[k-1,(EE-E0)/htot]*ChebyshevT[q-1,(EE-E0)/htot]*EE];
			Smat[[k,q]]=Total[Pro\[Psi]*ChebyshevT[k-1,(EE-E0)/htot]*ChebyshevT[q-1,(EE-E0)/htot]];
		),{q,1,d}]
	),{k,1,d}];
	Hmat=(Hmat+ConjugateTranspose[Hmat])/2.;
	Smat=(Smat+ConjugateTranspose[Smat])/2.;
	Return[{Hmat,Smat}]
)]


funMatGaussianPower[EE_,Pro\[Psi]_,d_,htot_,\[Tau]_,E0_]:=Module[{costList,ProG\[Psi],Hmat,Smat},(
	costList=funCost[htot,\[Tau],d];
	ProG\[Psi]=Exp[-(EE-E0)^2*\[Tau]^2]*Pro\[Psi];
	Hmat=Table[Table[0.,d],d];
	Smat=Table[Table[0.,d],d];
	Do[(
		Do[(
			If[k+q-2==0,(
				Hmat[[k,q]]=Total[ProG\[Psi]*EE]/(costList[[k]]*costList[[q]]);
				Smat[[k,q]]=Total[ProG\[Psi]]/(costList[[k]]*costList[[q]])
			),(
				Hmat[[k,q]]=Total[ProG\[Psi]*(EE-E0)^(k+q-2)*EE]/(costList[[k]]*costList[[q]]);
				Smat[[k,q]]=Total[ProG\[Psi]*(EE-E0)^(k+q-2)]/(costList[[k]]*costList[[q]])
			)]
		),{q,1,d}]
	),{k,1,d}];
	Hmat=(Hmat+ConjugateTranspose[Hmat])/2.;
	Smat=(Smat+ConjugateTranspose[Smat])/2.;
	Return[{Hmat,Smat}]
)]


(*funMatGaussianPower[EE_,Pro\[Psi]_,d_,\[Tau]_,E0_]:=Module[{ProG\[Psi],Hmat,Smat},(
	ProG\[Psi]=Exp[-(EE-E0)^2*\[Tau]^2]*Pro\[Psi];
	Hmat=Table[Table[0.,d],d];
	Smat=Table[Table[0.,d],d];
	Do[(
		Do[(
			If[k+q-2==0,(
				Hmat[[k,q]]=Total[ProG\[Psi]*EE]/(funCost[\[Tau],k]*funCost[\[Tau],q]);
				Smat[[k,q]]=Total[ProG\[Psi]]/(funCost[\[Tau],k]*funCost[\[Tau],q])
			),(
				Hmat[[k,q]]=Total[ProG\[Psi]*(EE-E0)^(k+q-2)*EE]/(funCost[\[Tau],k]*funCost[\[Tau],q]);
				Smat[[k,q]]=Total[ProG\[Psi]*(EE-E0)^(k+q-2)]/(funCost[\[Tau],k]*funCost[\[Tau],q])
			)]
		),{q,1,d}]
	),{k,1,d}];
	Hmat=(Hmat+ConjugateTranspose[Hmat])/2.;
	Smat=(Smat+ConjugateTranspose[Smat])/2.;
	Return[{Hmat,Smat}]
)]*)


funMatInversePower[EE_,Pro\[Psi]_,d_,E0_]:=Module[{Hmat,Smat},(
	Hmat=Table[Table[0.,d],d];
	Smat=Table[Table[0.,d],d];
	Do[(
		Do[(
			If[k+q-2==0,(
				Hmat[[k,q]]=Total[Pro\[Psi]*EE];
				Smat[[k,q]]=Total[Pro\[Psi]]
			),(
				Hmat[[k,q]]=Total[Pro\[Psi]*(EE-E0)^(-k-q+2)*EE];
				Smat[[k,q]]=Total[Pro\[Psi]*(EE-E0)^(-k-q+2)]
			)]
			
		),{q,1,d}]
	),{k,1,d}];
	Hmat=(Hmat+ConjugateTranspose[Hmat])/2.;
	Smat=(Smat+ConjugateTranspose[Smat])/2.;
	Return[{Hmat,Smat}]
)]


funMatITE[EE_,Pro\[Psi]_,d_,\[Tau]_,E0_]:=Module[{ITE,Hmat,Smat},(
	ITE=Exp[-(EE-E0)*\[Tau]];
	Hmat=Table[Table[0.,d],d];
	Smat=Table[Table[0.,d],d];
	Do[(
		Do[(
			If[k+q-2==0,(
				Hmat[[k,q]]=Total[Pro\[Psi]*EE];
				Smat[[k,q]]=Total[Pro\[Psi]]
			),(
				Hmat[[k,q]]=Total[Pro\[Psi]*ITE^(k+q-2)*EE];
				Smat[[k,q]]=Total[Pro\[Psi]*ITE^(k+q-2)]
			)]
		),{q,1,d}]
	),{k,1,d}];
	Hmat=(Hmat+ConjugateTranspose[Hmat])/2.;
	Smat=(Smat+ConjugateTranspose[Smat])/2.;
	Return[{Hmat,Smat}]
)]


funMatRTE[EE_,Pro\[Psi]_,d_,\[CapitalDelta]t_,E0_]:=Module[{RTE,Hmat,Smat},(
	RTE=Exp[I*(EE-E0)*\[CapitalDelta]t];
	Hmat=Table[Table[0.,d],d];
	Smat=Table[Table[0.,d],d];
	Do[(
		Do[(
			If[k+q-2==0,(
				Hmat[[k,q]]=Total[Pro\[Psi]*EE];
				Smat[[k,q]]=Total[Pro\[Psi]]
			),(
				Hmat[[k,q]]=Total[Pro\[Psi]*RTE^(k-q)*EE];
				Smat[[k,q]]=Total[Pro\[Psi]*RTE^(k-q)]
			)]
		),{q,1,d}]
	),{k,1,d}];
	Hmat=(Hmat+ConjugateTranspose[Hmat])/2.;
	Smat=(Smat+ConjugateTranspose[Smat])/2.;
	Return[{Hmat,Smat}]
)]


funMatFilter[EE_,Pro\[Psi]_,d_,T_,E0_,\[CapitalDelta]E_]:=Module[{ProG\[Psi],Hmat,Smat},(
	Hmat=Table[Table[0.,d],d];
	Smat=Table[Table[0.,d],d];
	Do[(
		Do[(
			ProF\[Psi]=Sinc[(EE-(E0+(k-1)*\[CapitalDelta]E))*T]*Sinc[(EE-(E0+(q-1)*\[CapitalDelta]E))*T]*Pro\[Psi];
			Hmat[[k,q]]=Total[ProF\[Psi]*EE];
			Smat[[k,q]]=Total[ProF\[Psi]]
		),{q,1,d}]
	),{k,1,d}];
	Hmat=(Hmat+ConjugateTranspose[Hmat])/2.;
	Smat=(Smat+ConjugateTranspose[Smat])/2.;
	Return[{Hmat,Smat}]
)]


(*Plot*)


funGammaEpsilon[Eg_,pg_,Hmat_,Smat_,Ide_,CH_,CS_,log\[Eta]List_]:=Module[{\[Gamma]List,\[Epsilon]List,log\[Eta],\[Eta],EK,cn,\[Epsilon],\[Gamma]},(
	\[Gamma]List=0.*log\[Eta]List;
	\[Epsilon]List=0.*log\[Eta]List;
	Do[(
		log\[Eta]=log\[Eta]List[[j]];
		\[Eta]=10.^log\[Eta];
		
		{EK,cn}=funDiagonalisation[Hmat+2.*CH*\[Eta]*Ide,Smat+2.*CS*\[Eta]*Ide];
		\[Epsilon]=EK-Eg;
		
		\[Gamma]=(pg^2\[Epsilon]^2)/(16\[Eta]^2);
		\[Gamma]List[[j]]=\[Gamma];
		\[Epsilon]List[[j]]=\[Epsilon];
		(*Print[{"j",j,\[Eta],\[Gamma],\[Epsilon],ToString[Now]}];*)
	),{j,1,Length[log\[Eta]List]}];
	Return[{\[Gamma]List,\[Epsilon]List}]
)]


funInterpolation[xList_,yList_,x_]:=Module[{i,y},(
	i=Position[(x-xList)^2,Min[(x-xList)^2]][[1,1]];
	If[(xList[[i]]<x&&i<Length[xList])||(xList[[i]]>x&&i==1),y=yList[[i]]+(x-xList[[i]]) (yList[[i+1]]-yList[[i]])/(xList[[i+1]]-xList[[i]])];
	If[xList[[i]]==x,y=yList[[i]]];
	If[(xList[[i]]>x&&i>1)||(xList[[i]]<x&&i==Length[xList]),y=yList[[i]]+(x-xList[[i]]) (yList[[i-1]]-yList[[i]])/(xList[[i-1]]-xList[[i]])];
	Return[y]
)]
