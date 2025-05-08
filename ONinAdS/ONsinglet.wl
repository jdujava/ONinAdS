(* ::Package:: *)

(* ::Section:: *)
(*Prolog*)


BeginPackage["ONinAdS`ONsinglet`", "ONinAdS`CFT`"];


(* clear all declared variables/functions which memoize their values *)
ClearAll[eq2d, eq4d];
ClearAll[poleSinglet2d, polesList2d, poleSinglet4d, polesList4d];
ClearAll[opeSqList2d, opeSqList4d];


(* expose these variables/functions to public *)
B; B2d; B4d; B6d; BubbleConformal;
Spec\[CapitalGamma]; Spec\[CapitalGamma]2d; Spec\[CapitalGamma]4d; SpecB; SpecB2d; SpecB4d; SpecExch; SpecExch2d; SpecExch4d;
SpecSinglet; SpecSinglet2d; SpecSinglet4d;
eq2d::usage = "eq2d[\[CapitalDelta]ext,\[Lambda]] returns the equation (in parameter \[CapitalDelta]) 2B+1/\[Lambda]==0 for d=2 and external scaling dimension \[CapitalDelta]ext";
eq4d::usage = "eq4d[\[CapitalDelta]ext,\[Lambda]] returns the equation (in parameter \[CapitalDelta]) 2B+1/\[Lambda]==0 for d=4 and external scaling dimension \[CapitalDelta]ext";
Nmax2d; Nmax4d; poleSinglet2d; poleSinglet4d;
polesList2d::usage = "polesList2d[\[CapitalDelta]ext,\[Lambda]] returns roots of 2B+1/\[Lambda]=0 (poles of SpecB)";
polesList4d::usage = "polesList4d[\[CapitalDelta]ext,\[Lambda]] returns roots of 2B+1/\[Lambda]=0 (poles of SpecB)";
opeSqList2d::usage = "opeSqList2d[nmax,\[CapitalDelta]ext,\[Lambda]] returns squared OPE coefficients for d=2";
opeSqList4d::usage = "opeSqList4d[nmax,\[CapitalDelta]ext,\[Lambda]] returns squared OPE coefficients for d=4";


Begin["`Private`"];


(* ::Section:: *)
(*Spectral function in the singlet sector*)


(* ::Subsection:: *)
(*Bubble function*)


(* Bubble function in general dimension (converges for d+1<4) *)
(* (4.27)[Carmi, di Pietro, Komatsu] *)
(* \[CapitalDelta]ext: external scaling dimension <-> mass of \[CapitalPhi] fields *)
(* \[CapitalDelta]=(d/2+I \[Nu]): spectral parameter *)
B[d_,\[CapitalDelta]ext_,\[CapitalDelta]_] := With[
	{s\[CapitalDelta]=d-\[CapitalDelta]},
	(Gamma[\[CapitalDelta]ext]Gamma[\[CapitalDelta]ext-d/2+1/2]Gamma[2\[CapitalDelta]ext-d/2])/(4(4\[Pi])^(d/2)) (
	  +Gamma[\[CapitalDelta]ext- \[CapitalDelta]/2] HypergeometricPFQRegularized[{d/2,\[CapitalDelta]ext,\[CapitalDelta]ext-d/2+1/2,\[CapitalDelta]ext- \[CapitalDelta]/2,2\[CapitalDelta]ext-d/2},{\[CapitalDelta]ext+1/2,\[CapitalDelta]ext-d/2+1,\[CapitalDelta]ext- \[CapitalDelta]/2+1,2\[CapitalDelta]ext-d+1},1]
	  +Gamma[\[CapitalDelta]ext-s\[CapitalDelta]/2] HypergeometricPFQRegularized[{d/2,\[CapitalDelta]ext,\[CapitalDelta]ext-d/2+1/2,\[CapitalDelta]ext-s\[CapitalDelta]/2,2\[CapitalDelta]ext-d/2},{\[CapitalDelta]ext+1/2,\[CapitalDelta]ext-d/2+1,\[CapitalDelta]ext-s\[CapitalDelta]/2+1,2\[CapitalDelta]ext-d+1},1]
	)
]


(* ::Subsubsection:: *)
(*Bootstrap calculation of the Bubble function*)


(* ::Code:: *)
(*(* Obtained by requiring cancelation of spurious J=0 MFT poles in the singlet sector *)*)
(*(* (3.26)[Carmi, di Pietro, Komatsu] *)*)
(*(* In Bsummand, \[CapitalDelta] corresponds to the dimension of external operators, and \[Nu] is the spectral parameter *)*)
(*(* Note, that in other places we usually use \[CapitalDelta]ext for the external dimension,*)
(*     and work directly with \[CapitalDelta]=(d/2+I \[Nu]) instead of \[Nu]*)*)
(*Bsummand = (2\[CapitalDelta]+2n-d/2)/(\[Nu]^2+(2\[CapitalDelta]+2n-d/2)^2) ( *)
(*	(Pochhammer[d/2,n]Gamma[\[CapitalDelta]+n]Gamma[\[CapitalDelta]+n-d/2+1/2]Gamma[2\[CapitalDelta]+n-d/2])*)
(*	/((4\[Pi])^(d/2) Gamma[n+1]Gamma[\[CapitalDelta]+n+1/2]Gamma[\[CapitalDelta]+n-d/2+1]Gamma[2\[CapitalDelta]-d+n+1])*)
(*);*)


(* ::Code:: *)
(*(* direct sum gives expression which Mathematica is unable to simplify further *)*)
(*Sum[FullSimplify[Bsummand], {n,0,\[Infinity]}] // FullSimplify*)
(*(* saved result from direct evaluation of the sum *)*)
(*BsumDirect = -(1/\[Nu])I 2^(-2-d) \[Pi]^(-d/2) Gamma[1-d/2+2\[CapitalDelta]] ( *)
(*	  d Gamma[1+\[CapitalDelta]] Gamma[3/2-d/2+\[CapitalDelta]] ( *)
(*		  Gamma[1-d/4+\[CapitalDelta]-(I \[Nu])/2] HypergeometricPFQRegularized[{1+d/2,1+\[CapitalDelta],3/2-d/2+\[CapitalDelta],1-d/2+2 \[CapitalDelta],1-d/4+\[CapitalDelta]-(I \[Nu])/2},{3/2+\[CapitalDelta],2-d/2+\[CapitalDelta],2-d+2\[CapitalDelta],2-d/4+\[CapitalDelta]-(I \[Nu])/2},1]*)
(*		- Gamma[1-d/4+\[CapitalDelta]+(I \[Nu])/2] HypergeometricPFQRegularized[{1+d/2,1+\[CapitalDelta],3/2-d/2+\[CapitalDelta],1-d/2+2 \[CapitalDelta],1-d/4+\[CapitalDelta]+(I \[Nu])/2},{3/2+\[CapitalDelta],2-d/2+\[CapitalDelta],2-d+2\[CapitalDelta],2-d/4+\[CapitalDelta]+(I \[Nu])/2},1]*)
(*		)*)
(*	+ Gamma[\[CapitalDelta]] Gamma[1/2-d/2+\[CapitalDelta]] ( *)
(*		  Gamma[-(d/4)+\[CapitalDelta]-(I \[Nu])/2] HypergeometricPFQRegularized[{d/2,\[CapitalDelta],1/2-d/2+\[CapitalDelta],-(d/2)+2\[CapitalDelta],-(d/4)+\[CapitalDelta]-(I \[Nu])/2},{1/2+\[CapitalDelta],1-d/2+\[CapitalDelta],1-d+2\[CapitalDelta],1-d/4+\[CapitalDelta]-(I \[Nu])/2},1]*)
(*		- Gamma[-(d/4)+\[CapitalDelta]+(I \[Nu])/2] HypergeometricPFQRegularized[{d/2,\[CapitalDelta],1/2-d/2+\[CapitalDelta],-(d/2)+2\[CapitalDelta],-(d/4)+\[CapitalDelta]+(I \[Nu])/2},{1/2+\[CapitalDelta],1-d/2+\[CapitalDelta],1-d+2\[CapitalDelta],1-d/4+\[CapitalDelta]+(I \[Nu])/2},1]*)
(*	)*)
(*);*)


(* ::Code:: *)
(*(* now we will proceed by calculating the part for \[Nu]->0, and the rest separately, summing only at the end *)*)
(*(* this is valid only for d+1<4, where the sum converges *)*)
(*B\[Nu]0 = Sum[FullSimplify[Bsummand/.\[Nu]->0], {n,0,\[Infinity]}]*)
(*BsumSubtr = FullSimplify[Sum[FullSimplify[Bsummand - (Bsummand /. \[Nu]->0)], {n,0,\[Infinity]}]]*)
(*BsumFull = FullSimplify[BsumSubtr + B\[Nu]0]*)


(* ::Code:: *)
(*(* it really is the same as the definition at the start of the subsection *)*)
(*B[d, \[CapitalDelta], d/2+I \[Nu]]-BsumFull // FullSimplify*)


(* ::Code:: *)
(*(* alternatively, we can just calculate with the first summand, and restore shadow symmetry at the end *)*)
(*Apart[Bsummand/.{\[CapitalDelta]->\[CapitalDelta]ext,\[Nu]->-I(-d/2+\[CapitalDelta])},\[CapitalDelta]]*)
(*Sum[FullSimplify[First[%]], {n,0,\[Infinity]}] // FullSimplify*)
(*(* it indeed gives the same result *)*)
(*(% + (% /. \[CapitalDelta]->d-\[CapitalDelta])) - B[d,\[CapitalDelta]ext,\[CapitalDelta]] // FullSimplify*)
(**)


(* ::Subsubsection::Closed:: *)
(*Specially for d=2*)


B2d[\[CapitalDelta]ext_,\[CapitalDelta]_] := With[
	{\[Nu]=(\[CapitalDelta]-2/2)/I},
	I/(8 \[Pi] \[Nu]) (PolyGamma[\[CapitalDelta]ext-(1+I \[Nu])/2] - PolyGamma[\[CapitalDelta]ext-(1-I \[Nu])/2])
]


B2dAlt[\[CapitalDelta]ext_,\[CapitalDelta]_] := With[
	{s\[CapitalDelta]=2-\[CapitalDelta]},
	 1/(4\[Pi])(PolyGamma[\[CapitalDelta]ext-\[CapitalDelta]/2] - PolyGamma[\[CapitalDelta]ext-s\[CapitalDelta]/2])/(s\[CapitalDelta]-\[CapitalDelta])
]


(* ::Code:: *)
(*(* coincides with d->2 substitution in general B *)*)
(*FullSimplify[B[2,\[CapitalDelta]ext,\[CapitalDelta]]]*)
(*B2d[\[CapitalDelta]ext,\[CapitalDelta]] - FullSimplify[B[2,\[CapitalDelta]ext,\[CapitalDelta]]] // FullSimplify*)
(*B2d[\[CapitalDelta]ext,\[CapitalDelta]] - B2dAlt[\[CapitalDelta]ext,\[CapitalDelta]] // FullSimplify*)


(* ::Code:: *)
(*Plot[B2d[1,\[CapitalDelta]],{\[CapitalDelta],1,20},PlotRange->{{1,20},{-0.3,0.3}}]*)


(* ::Code:: *)
(*(* value at \[CapitalDelta]=d/2=1 *)*)
(*Limit[B2d[\[CapitalDelta]ext,\[CapitalDelta]],\[CapitalDelta]->1]*)
(*%/.\[CapitalDelta]ext->1*)


(* ::Code:: *)
(*(* for integer \[CapitalDelta]ext, the expression further simplifies *)*)
(*Table[B2d[\[CapitalDelta]ext,\[CapitalDelta]] // FunctionExpand // FullSimplify, {\[CapitalDelta]ext,1,4}]*)


(* ::Subsubsection::Closed:: *)
(*Specially for d=4 *)


(* ::Code:: *)
(*(* to obtain nice expression for d=4, we must subtitute the dimension*)
(*	already before calculating the sums *)*)
(*(* NOTE: in d=4 there is a divergence, regularization is important *)*)
(*Bsummand4d = FullSimplify[Bsummand/.d->4];*)
(*BsumSubtr4d = FullSimplify[Sum[Simplify[Bsummand4d -(Bsummand4d /. \[Nu]->0)], {n,0,\[Infinity]}]];*)
(*(* we fix the subtraction ambiguity by requirement of conformality at (\[CapitalDelta]ext=d-1=3, \[Lambda]=\[Infinity]), as shown below *)*)
(*BsumFull4dFixedbyConf = BsumSubtr4d - 1/(16 \[Pi]^2)*)


B4d[\[CapitalDelta]ext_,\[CapitalDelta]_] := With[
	{\[Nu]=(\[CapitalDelta]-4/2)/I},
	-1/(16 \[Pi]^2)+
	\[Nu]/(128 \[Pi]^2 (1+\[Nu]^2)) (
	  2 (-5+2 \[CapitalDelta]ext) \[Nu]
	  - I (4(-2+\[CapitalDelta]ext)^2 + \[Nu]^2) FunctionExpand@(HarmonicNumber[-2+\[CapitalDelta]ext-(I \[Nu])/2] - HarmonicNumber[-2+\[CapitalDelta]ext+(I \[Nu])/2])
	)
]


(* ::Code:: *)
(*Plot[B4d[3,\[CapitalDelta]],{\[CapitalDelta],2,20}]*)


(* ::Code:: *)
(*(* value at \[CapitalDelta]=d/2=2 -> thats how we chose to fix the subtraction ambiguity *)*)
(*Limit[B4d[\[CapitalDelta]ext,\[CapitalDelta]],\[CapitalDelta]->2]*)


(* ::Code:: *)
(*(* for integer \[CapitalDelta]ext, the expression further simplifies *)*)
(*Table[B4d[\[CapitalDelta]ext,\[CapitalDelta]] // FullSimplify, {\[CapitalDelta]ext,2,5}]*)


(* ::Subsubsection::Closed:: *)
(*Large \[CapitalDelta] expansions*)


(* ::Code:: *)
(*$Assumptions={\[CapitalDelta]\[Element]Reals,\[CapitalDelta]>0}*)


(* ::Code:: *)
(*Series[B2d[\[CapitalDelta]ext,\[CapitalDelta]],{\[CapitalDelta],\[Infinity],2}] // Normal // FullSimplify*)


(* ::Code:: *)
(*Series[B4d[\[CapitalDelta]ext,\[CapitalDelta]],{\[CapitalDelta],\[Infinity],0}] // Normal // FullSimplify*)
(*Series[B4d[\[CapitalDelta]ext,\[CapitalDelta]]/\[CapitalDelta],{\[CapitalDelta],\[Infinity],0}] // Normal // FullSimplify*)


(* ::Subsubsection::Closed:: *)
(*Specially for d=6*)


(* ::Code:: *)
(*(* to obtain nice expression for d=6, we must subtitute the dimension*)
(*	already before calculating the sums *)*)
(*(* NOTE: in d=6 there is a divergence, regularization is important *)*)
(*Bsummand6d = FullSimplify[Bsummand/.d->6]*)
(*subtraction=Series[Bsummand6d,{\[Nu],0,2}]//Normal*)
(*FullSimplify[Bsummand6d - subtraction]*)
(*BsumSubtr6d = FullSimplify[Sum[%, {n,0,\[Infinity]}]]*)


(* saved result of BsumSubtr6d, with added piece such that \[CapitalDelta]ext=d-1=5 is conformal (as shown in the following) *)
B6d[\[CapitalDelta]ext_,\[CapitalDelta]_] := With[
	{\[Nu]=(\[CapitalDelta]-6/2)/I},
	1/(393216 \[Pi]^3 \[Nu] (4+\[Nu]^2)) (24 \[Nu] (16 (7-2 \[CapitalDelta]ext)^2+32 (25+2 (-7+\[CapitalDelta]ext) \[CapitalDelta]ext) \[Nu]^2+(-227+3 \[CapitalDelta]ext (99+4 (-9+\[CapitalDelta]ext) \[CapitalDelta]ext)) \[Nu]^4)+96 (-5+2 \[CapitalDelta]ext+I \[Nu]) (-I+\[Nu]) (I+\[Nu]) (-5 I+2 I \[CapitalDelta]ext+\[Nu]) ((7-2 \[CapitalDelta]ext)^2+\[Nu]^2) PolyGamma[0,-(3/2)+\[CapitalDelta]ext-(I \[Nu])/2]-96 I (1+\[Nu]^2) ((5-2 \[CapitalDelta]ext)^2+\[Nu]^2) ((7-2 \[CapitalDelta]ext)^2+\[Nu]^2) PolyGamma[0,-(3/2)+\[CapitalDelta]ext+(I \[Nu])/2]+\[Nu] (4+\[Nu]^2) (-6 (4 (35-24 \[CapitalDelta]ext+4 \[CapitalDelta]ext^2)^2+(3971+8 (-6+\[CapitalDelta]ext) \[CapitalDelta]ext (109+6 (-6+\[CapitalDelta]ext) \[CapitalDelta]ext)) \[Nu]^2) PolyGamma[1,-(5/2)+\[CapitalDelta]ext]+(35-24 \[CapitalDelta]ext+4 \[CapitalDelta]ext^2)^2 \[Nu]^2 PolyGamma[3,-(5/2)+\[CapitalDelta]ext]))
	(*-(\[Pi] (6399-75 \[Pi]^2 (-3+\[CapitalDelta])^2+811 (-6+\[CapitalDelta]) \[CapitalDelta]))/(131072 \[Pi]^2)*) (* extra piece *)
	+(900 - (-811 + 75 \[Pi]^2) \[Nu]^2)/(131072 \[Pi]) (* extra piece *)
]


(* ::Code:: *)
(*Plot[B6d[5,\[CapitalDelta]],{\[CapitalDelta],3,20}]*)


(* ::Subsection::Closed:: *)
(*Conformal points*)


(* ::Code:: *)
(*(* let's evaluate bubble functions for \[CapitalDelta]ext=d-1, which should correspond to critical points (for \[Lambda]=\[Infinity]) *)*)
(*(* we obtain very special forms (with appropriately chosen subtractions) *)*)
(*B2dConf = B2d[1,\[CapitalDelta]] // FullSimplify*)
(*B4dConf = B4d[3,\[CapitalDelta]] // FullSimplify*)
(*B6dConf = B6d[5,\[CapitalDelta]] // FullSimplify*)


(* ::Code:: *)
(*(* try to guess general formula for even d *)*)
(*Ratios[{8,128,4096,196608}]*)
(*Table[8 16^(d/2-1)(d/2-1)!,{d,2,8,2}]*)
(*8 16^(d/2-1)(d/2-1)! //FullSimplify*)


(* below, we will check this also for d=8, d=10, and d=12 *)
BubbleConformal[d_,\[CapitalDelta]_] := - 1/(2^(-1+2 d) \[Pi]^(d/2-1) Gamma[d/2]) (Product[\[CapitalDelta]-i,{i,4-d,2(d-2),2}]/Product[\[CapitalDelta]-i,{i,1,d-1,2}]) Cot[\[Pi]/2 \[CapitalDelta]]
BubbleConformalAlt[d_,\[CapitalDelta]_] := - 1/(2^(d+2) \[Pi]^(d/2-1) Gamma[d/2])(Gamma[(\[CapitalDelta]-(4-d))/2+1]/Gamma[(\[CapitalDelta]-2(d-2))/2])/(Gamma[(\[CapitalDelta]+1)/2]/Gamma[(-(d-\[CapitalDelta])+1)/2]) Cot[\[Pi]/2 \[CapitalDelta]]


(* ::Code:: *)
(*(* check that the formulas indeed reproduce d=2,d=4,d=6 results *)*)
(*Table[BubbleConformal[d,\[CapitalDelta]],{d,2,10,2}]*)
(*Table[BubbleConformalAlt[d,\[CapitalDelta]],{d,2,10,2}] // FunctionExpand;*)
(*%%-%*)


(* ::Code:: *)
(*Plot[B2dConf,{\[CapitalDelta],1,16},PlotRange->{-1/4,1/4}]*)
(*Plot[B4dConf,{\[CapitalDelta],2,16}]*)
(*Plot[B6dConf,{\[CapitalDelta],3,16}]*)


(* ::Code:: *)
(*(* the values at integer \[CapitalDelta] show special behavior *)*)
(*Limit[{#,B2dConf},\[CapitalDelta]->#]&/@Range[1,16]*)
(*Limit[{#,B4dConf},\[CapitalDelta]->#]&/@Range[2,16]*)
(*Limit[{#,B6dConf},\[CapitalDelta]->#]&/@Range[3,16] // FullSimplify*)


(* ::Subsubsection::Closed:: *)
(*Conformal point for d=8*)


(* ::Code:: *)
(*(* need to define Bsummand (in initial bootstrap calculation of bubble function) *)*)
(*(* NOTE: in d=8 there is a divergence, regularization is important (up to 4th order in \[Nu]) *)*)
(*Bsummand8d = FullSimplify[Bsummand/.d->8]*)
(*subtraction=Series[Bsummand8d,{\[Nu],0,4}]//Normal*)
(*BsumSubtr8d = Sum[Bsummand8d - subtraction, {n,0,\[Infinity]}]*)


(* saved result of BsumSubtr8d *)
BsumSubtr8d[\[CapitalDelta]ext_,\[CapitalDelta]_] := With[
	{\[Nu]=(\[CapitalDelta]-8/2)/I},
	-(((-3+\[CapitalDelta]ext)^2 (-1+\[CapitalDelta]ext) \[Nu]^6 HypergeometricPFQ[{4,-(7/2)+\[CapitalDelta]ext,-2+\[CapitalDelta]ext,-2+\[CapitalDelta]ext,-2+\[CapitalDelta]ext,-2+\[CapitalDelta]ext,-2+\[CapitalDelta]ext,\[CapitalDelta]ext,-4+2 \[CapitalDelta]ext,-2+\[CapitalDelta]ext-(I \[Nu])/2,-2+\[CapitalDelta]ext+(I \[Nu])/2},{-3+\[CapitalDelta]ext,-1+\[CapitalDelta]ext,-1+\[CapitalDelta]ext,-1+\[CapitalDelta]ext,-1+\[CapitalDelta]ext,-1+\[CapitalDelta]ext,1/2+\[CapitalDelta]ext,-7+2 \[CapitalDelta]ext,-1+\[CapitalDelta]ext-(I \[Nu])/2,-1+\[CapitalDelta]ext+(I \[Nu])/2},1])/(256 \[Pi]^4 (-2+\[CapitalDelta]ext)^4 (3+4 (-2+\[CapitalDelta]ext) \[CapitalDelta]ext) (4 (-2+\[CapitalDelta]ext)^2+\[Nu]^2)))
]


(* ::Code:: *)
(*(* since Mathematica is unable to effeciently work with B8d[7,\[CapitalDelta]], calculate it directly*)
(*   by setting \[CapitalDelta]ext=d-1=7 at the level of summand, as it is then able to get simplified result *)*)
(*Bsummand8dconf = FullSimplify[Bsummand/.{d->8,\[CapitalDelta]->7}]*)
(*subtraction=Series[Bsummand8dconf,{\[Nu],0,4}]//Normal;*)
(*BsummandSubtracted=Bsummand8dconf - subtraction // FullSimplify*)
(*shouldbezeros=ParallelTable[Sum[Simplify@(BsummandSubtracted/.\[Nu]->(\[CapitalDelta]-4)/I), {n,0,\[Infinity]}],{\[CapitalDelta],9,13,2}]*)


(* ::Code:: *)
(*(* fix subtractions by requiring zeros at odd integer positions *)
(*     (corresponding to the displacement operator and after)  *)*)
(*finitepiece8d=a-b (\[CapitalDelta]-4)^2+c (\[CapitalDelta]-4)^4;*)
(*equations8d=(#==0)&/@(Table[finitepiece8d,{\[CapitalDelta],9,13,2}]+shouldbezeros)*)


(* ::Code:: *)
(*finitepiecesolution=First@Solve[equations8d,{a,b,c}]*)
(*finpiece8d=finitepiece8d/.finitepiecesolution // Simplify*)


(* ::Code:: *)
(*(* check that we indeed obtain zeros even after the first ones fixed by subtractions *)*)
(*ParallelTable[finpiece8d+Sum[Simplify@(BsummandSubtracted/.\[Nu]->(\[CapitalDelta]-4)/I), {n,0,\[Infinity]}],{\[CapitalDelta],9,23,2}] // Simplify*)


(* ::Code:: *)
(*(* values at integer positions can be explicitly evaluated *)*)
(*ParallelTable[finpiece8d+Sum[Simplify@(BsummandSubtracted/.\[Nu]->(\[CapitalDelta]-4)/I), {n,0,\[Infinity]}],{\[CapitalDelta],4,10}] // Simplify*)


(* ::Code:: *)
(*(* comparing with the "general" formula, we obtain precise match *)*)
(*guess8dconf=BubbleConformal[8,\[CapitalDelta]]*)
(*Limit[guess8dconf,\[CapitalDelta]->#]&/@Range[4,10]*)
(*Plot[guess8dconf,{\[CapitalDelta],4,20},PlotRange->{-0.2,0.2}]*)


(* ::Code:: *)
(*(* numerical check for noninteger \[CapitalDelta] *)*)
(*\[CapitalDelta]sub=5.2;*)
(*BubblebySum=finpiece8d+ParallelSum[Simplify@(BsummandSubtracted/.\[Nu]->(\[CapitalDelta]sub-4)/I), {n,0,10^3}] /.\[CapitalDelta]->\[CapitalDelta]sub //EchoTiming*)
(*BubblebyGuess=guess8dconf /.\[CapitalDelta]->\[CapitalDelta]sub*)
(*BubblebySum-BubblebyGuess*)


(* ::Subsubsection::Closed:: *)
(*Conformal point for d=10*)


(* ::Code:: *)
(*(* need to define Bsummand (in initial bootstrap calculation of bubble function) *)*)
(*(* NOTE: in d=10 there is a divergence, regularization is important (up to 6th order in \[Nu]) *)*)
(*Bsummand10dconf = FullSimplify[Bsummand/.{d->10,\[CapitalDelta]->9}]*)
(*subtraction=Series[Bsummand10dconf,{\[Nu],0,6}]//Normal;*)
(*BsummandSubtracted=Bsummand10dconf - subtraction // FullSimplify*)
(*shouldbezeros=ParallelTable[Sum[Simplify@(BsummandSubtracted/.\[Nu]->(\[CapitalDelta]-5)/I), {n,0,\[Infinity]}],{\[CapitalDelta],11,17,2}]*)


(* ::Code:: *)
(*(* fix subtractions by requiring zeros at odd integer positions *)
(*     (corresponding to the displacement operator and after)  *)*)
(*finitepiece10d=a-b (\[CapitalDelta]-5)^2+c (\[CapitalDelta]-5)^4-d (\[CapitalDelta]-5)^6;*)
(*equations10d=(#==0)&/@(Table[finitepiece10d,{\[CapitalDelta],11,17,2}]+shouldbezeros)*)


(* ::Code:: *)
(*finitepiecesolution=First@Solve[equations10d,{a,b,c,d}]*)
(*finpiece10d=finitepiece10d/.finitepiecesolution // Simplify*)


(* ::Code:: *)
(*(* check that we indeed obtain zeros even after the first ones fixed by subtractions *)*)
(*ParallelTable[finpiece10d+Sum[Simplify@(BsummandSubtracted/.\[Nu]->(\[CapitalDelta]-5)/I), {n,0,\[Infinity]}],{\[CapitalDelta],11,21,2}] // Simplify*)


(* ::Code:: *)
(*(* values at integer positions can be explicitly evaluated *)*)
(*ParallelTable[finpiece10d+Sum[Simplify@(BsummandSubtracted/.\[Nu]->(\[CapitalDelta]-5)/I), {n,0,\[Infinity]}],{\[CapitalDelta],5,13}] // Simplify*)


(* ::Code:: *)
(*(* comparing with the "general" formula, we obtain precise match *)*)
(*guess10dconf=BubbleConformal[10,\[CapitalDelta]]*)
(*Limit[guess10dconf,\[CapitalDelta]->#]&/@Range[5,13]*)
(*Plot[guess10dconf,{\[CapitalDelta],5,24},PlotRange->{-0.003,0.003}]*)
(*Plot[guess10dconf,{\[CapitalDelta],5,24},PlotRange->{-0.2,0.2}]*)


(* ::Code:: *)
(*(* numerical check for noninteger \[CapitalDelta] *)*)
(*\[CapitalDelta]sub=13.78;*)
(*BubblebySum=finpiece10d+ParallelSum[Simplify@(BsummandSubtracted/.\[Nu]->(\[CapitalDelta]sub-5)/I), {n,0,10^3}] /.\[CapitalDelta]->\[CapitalDelta]sub //EchoTiming*)
(*BubblebyGuess=guess10dconf /.\[CapitalDelta]->\[CapitalDelta]sub*)
(*BubblebySum-BubblebyGuess*)


(* ::Subsubsection::Closed:: *)
(*Conformal point for d=12*)


(* ::Code:: *)
(*(* need to define Bsummand (in initial bootstrap calculation of bubble function) *)*)
(*(* NOTE: in d=12 there is a divergence, regularization is important (up to 8th order in \[Nu]) *)*)
(*Bsummand12dconf = FullSimplify[Bsummand/.{d->12,\[CapitalDelta]->11}]*)
(*subtraction=Series[Bsummand12dconf,{\[Nu],0,8}]//Normal;*)
(*BsummandSubtracted=Bsummand12dconf - subtraction // FullSimplify*)
(*shouldbezeros=ParallelTable[Sum[Simplify@(BsummandSubtracted/.\[Nu]->(\[CapitalDelta]-6)/I), {n,0,\[Infinity]}],{\[CapitalDelta],13,21,2}]*)


(* ::Code:: *)
(*(* fix subtractions by requiring zeros at odd integer positions *)
(*     (corresponding to the displacement operator and after)  *)*)
(*finitepiece12d=a-b (\[CapitalDelta]-6)^2+c (\[CapitalDelta]-6)^4-d (\[CapitalDelta]-6)^6+e (\[CapitalDelta]-6)^8;*)
(*equations12d=(#==0)&/@(Table[finitepiece12d,{\[CapitalDelta],13,21,2}]+shouldbezeros);*)


(* ::Code:: *)
(*finitepiecesolution=First@Solve[equations12d,{a,b,c,d,e}]*)
(*finpiece12d=finitepiece12d/.finitepiecesolution // Simplify;*)


(* ::Code:: *)
(*(* check that we indeed obtain zeros even after the first ones fixed by subtractions *)*)
(*ParallelTable[finpiece12d+Sum[Simplify@(BsummandSubtracted/.\[Nu]->(\[CapitalDelta]-6)/I), {n,0,\[Infinity]}],{\[CapitalDelta],13,23,2}] // Simplify*)


(* ::Code:: *)
(*(* values at integer positions can be explicitly evaluated *)*)
(*ParallelTable[finpiece12d+Sum[Simplify@(BsummandSubtracted/.\[Nu]->(\[CapitalDelta]-6)/I), {n,0,\[Infinity]}],{\[CapitalDelta],6,15}] // Simplify*)


(* ::Code:: *)
(*(* comparing with the "general" formula, we obtain precise match *)*)
(*guess12dconf=BubbleConformal[12,\[CapitalDelta]]*)
(*Limit[guess12dconf,\[CapitalDelta]->#]&/@Range[6,15]*)
(*Plot[guess12dconf,{\[CapitalDelta],6,26},PlotRange->{-0.005,0.005}]*)
(*Plot[guess12dconf,{\[CapitalDelta],6,26},PlotRange->{-0.2,0.2}]*)


(* ::Code:: *)
(*(* numerical check for noninteger \[CapitalDelta] *)*)
(*\[CapitalDelta]sub=20.78;*)
(*BubblebySum=finpiece12d+ParallelSum[Simplify@(BsummandSubtracted/.\[Nu]->(\[CapitalDelta]sub-6)/I), {n,0,10^3}] /.\[CapitalDelta]->\[CapitalDelta]sub //EchoTiming*)
(*BubblebyGuess=guess12dconf /.\[CapitalDelta]->\[CapitalDelta]sub*)
(*BubblebySum-BubblebyGuess*)


(* ::Subsubsection::Closed:: *)
(*Odd dimensions*)


(* ::Code:: *)
(*(* for odd dimensions we do not see such behavior *)*)
(*(* in particular, the wouldbe operator following displacement operatr*)
(*   with dimensions \[CapitalDelta]=d+1+2n are in conflict with appearance of MFT poles *)
(*   at \[CapitalDelta]=2\[CapitalDelta]ext+2n, also at even \[CapitalDelta] *)*)


(* ::Text:: *)
(*Specially for d=3*)


(* ::Code:: *)
(*(* NOTE: in d=3 there is a divergence, regularization is important *)*)
(*Bsummand3d = FullSimplify[(Bsummand/.\[CapitalDelta]->d-1)/.d->3]*)
(*subtraction=Bsummand3d/.\[Nu]->0*)
(*FullSimplify[Bsummand3d - subtraction]*)
(*BsumSubtr3d = FullSimplify[Sum[%, {n,0,\[Infinity]}]/.\[Nu]->(\[CapitalDelta]-3/2)/I]*)


(* ::Code:: *)
(*B3dConf=BsumSubtr3d // FullSimplify*)
(*Plot[%,{\[CapitalDelta],3/2,20}]*)


(* ::Text:: *)
(*Specially for d=5*)


(* ::Code:: *)
(*(* to obtain nice expression for d=5, we must subtitute the dimension*)
(*	already before calculating the sums *)*)
(*(* NOTE: in d=5 there is a divergence, regularization is important *)*)
(*Bsummand5d = FullSimplify[(Bsummand/.\[CapitalDelta]->d-1)/.d->5]*)
(*subtraction=Series[Bsummand5d,{\[Nu],0,2}]//Normal*)
(*FullSimplify[Bsummand5d - subtraction]*)
(*BsumSubtr5d = FullSimplify[Sum[%, {n,0,\[Infinity]}]/.\[Nu]->(\[CapitalDelta]-5/2)/I]*)


(* ::Code:: *)
(*Limit[{#,BsumSubtr5d},\[CapitalDelta]->#]&/@Range[2,16] // FullSimplify*)


(* ::Code:: *)
(*B5dConf=BsumSubtr5d // FullSimplify*)
(*Plot[%,{\[CapitalDelta],3/2,20}]*)


(* ::Subsubsection::Closed:: *)
(*Noninteger d*)


(* ::Code:: *)
(*(* even for general dimension 2<d+1<4 and \[CapitalDelta]ext=d-1, we find first zero*)
(*	at the value corresponding to the displacement operator with \[CapitalDelta]=d+1 *)*)
(*Bconf = B[d,d-1,\[CapitalDelta]] // FullSimplify*)
(*ParallelTable[{\[CapitalDelta],FullSimplify@Bconf},{d,1.2,1.8,0.2},{\[CapitalDelta],d/2,6,0.1}];*)
(*ListPlot[#]&/@%*)


(* ::Code:: *)
(*(* the "general" formula probably does not really make sense for noneven d *)*)
(*Table[Plot[BubbleConformalAlt[d,\[CapitalDelta]],{\[CapitalDelta],d/2,10}] // FunctionExpand,{d,3,4,0.1}]*)


(* ::Subsection::Closed:: *)
(*Spectral functions *)


(* ::Subsubsection::Closed:: *)
(*Spectral function for J=0 s-channel \[Sigma]-exchange*)


(* from (4.18) [Carmi;di Pietro;Komatsu] *)
(* Note: we divide by K[s\[CapitalDelta]], since (4.18) is with respect to CBs,
         but our spectral function is with respect to CPWs *)
SpecExch[d_,\[CapitalDelta]ext_,\[Lambda]_,\[CapitalDelta]_] := Spec\[CapitalGamma][d,\[CapitalDelta]ext,\[CapitalDelta]] SpecB[d,\[CapitalDelta]ext,\[Lambda],\[CapitalDelta]]/K[d,d-\[CapitalDelta],0];
(* part including Gamma and other factors *)
Spec\[CapitalGamma][d_,\[CapitalDelta]ext_,\[CapitalDelta]_] := With[
	{s\[CapitalDelta]=d-\[CapitalDelta]},
	(Gamma[\[CapitalDelta]ext-\[CapitalDelta]/2]^2 Gamma[\[CapitalDelta]ext-s\[CapitalDelta]/2]^2 Gamma[\[CapitalDelta]/2]^4)
	/(4\[Pi]^(d/2) Gamma[\[CapitalDelta]ext]^2 Gamma[1-d/2+\[CapitalDelta]ext]^2 Gamma[\[CapitalDelta]-d/2] Gamma[\[CapitalDelta]])
]
(* part including the Bubble function and coupling \[Lambda] *)
SpecB[d_,\[CapitalDelta]ext_,\[Lambda]_,\[CapitalDelta]_] := -1/(1/\[Lambda]+2B[d,\[CapitalDelta]ext,\[CapitalDelta]]);


(* ::Code:: *)
(*Spec\[CapitalGamma][d,\[CapitalDelta]ext,\[CapitalDelta]]/K[d,d-\[CapitalDelta],0] // Simplify*)


(* ::Text:: *)
(*Specially for d=2*)


SpecB2d[\[CapitalDelta]ext_,\[Lambda]_,\[CapitalDelta]_] := -1/(1/\[Lambda] + 2 B2d[\[CapitalDelta]ext,\[CapitalDelta]]);
Spec\[CapitalGamma]2d[\[CapitalDelta]ext_,\[CapitalDelta]_] := Spec\[CapitalGamma][2,\[CapitalDelta]ext,\[CapitalDelta]];
SpecExch2d[\[CapitalDelta]ext_,\[Lambda]_,\[CapitalDelta]_] := SpecB2d[\[CapitalDelta]ext,\[Lambda],\[CapitalDelta]] Spec\[CapitalGamma]2d[\[CapitalDelta]ext,\[CapitalDelta]]/K[2,2-\[CapitalDelta],0];


(* ::Text:: *)
(*Specially for d=4 *)


SpecExch4d[\[CapitalDelta]ext_,\[Lambda]_,\[CapitalDelta]_] := SpecB4d[\[CapitalDelta]ext,\[Lambda],\[CapitalDelta]] Spec\[CapitalGamma]4d[\[CapitalDelta]ext,\[CapitalDelta]]/K[4,4-\[CapitalDelta],0];
SpecB4d[\[CapitalDelta]ext_,\[Lambda]_,\[CapitalDelta]_] := -1/(1/\[Lambda] + 2 B4d[\[CapitalDelta]ext,\[CapitalDelta]]);
Spec\[CapitalGamma]4d[\[CapitalDelta]ext_,\[CapitalDelta]_] := Spec\[CapitalGamma][4,\[CapitalDelta]ext,\[CapitalDelta]];


(* ::Subsubsection::Closed:: *)
(*Full O(1) piece of singlet spectral function \[LongDash] sum of t- and u-disconnected diagrams with the s-channel exchange*)


SpecSinglet[d_,\[CapitalDelta]ext_,\[Lambda]_,\[CapitalDelta]_,J_] := (1+(-1)^J)SpecDisc[d,\[CapitalDelta]ext,\[CapitalDelta],J] + KroneckerDelta[J,0]SpecExch[d,\[CapitalDelta]ext,\[Lambda],\[CapitalDelta]];


(* ::Text:: *)
(*d=2*)


SpecSinglet2d[\[CapitalDelta]ext_,\[Lambda]_,\[CapitalDelta]_,J_] := Simplify[(1+(-1)^J)SpecDisc[2,\[CapitalDelta]ext,\[CapitalDelta],J]] + KroneckerDelta[J,0]SpecExch2d[\[CapitalDelta]ext,\[Lambda],\[CapitalDelta]]


(* ::Code:: *)
(*SpecSinglet[2,\[CapitalDelta]ext,\[Lambda],\[CapitalDelta],J]-SpecSinglet2d[\[CapitalDelta]ext,\[Lambda],\[CapitalDelta],J] // FullSimplify*)


(* ::Text:: *)
(*d=4*)


SpecSinglet4d[\[CapitalDelta]ext_,\[Lambda]_,\[CapitalDelta]_,J_] := Simplify[(1+(-1)^J)SpecDisc[4,\[CapitalDelta]ext,\[CapitalDelta],J]] + KroneckerDelta[J,0]SpecExch4d[\[CapitalDelta]ext,\[Lambda],\[CapitalDelta]]


(* ::Section:: *)
(*Extracting data from the singlet spectral function*)


(* ::Subsection:: *)
(*simple poles of singlet spectral function \[DoubleLeftRightArrow] scaling dimensions of singlet operators*)


(* ::Code:: *)
(*(* simple poles = intersection points of the following graphs *)*)
(*Plot[*)
(*	{Evaluate[FullSimplify[2 B2d[1,1+z]]], -1/\[Lambda] /. \[Lambda]->5}, {z,0,30},*)
(*	ImageSize->700, AspectRatio->0.3, MaxRecursion->10, ExclusionsStyle->Directive[Gray,Opacity[0.4]]*)
(*]*)
(*Plot[*)
(*	{Evaluate[FullSimplify[2 B4d[3,2+z]]], -1/\[Lambda] /. \[Lambda]->31}, {z,0,30},*)
(*	ImageSize->700, AspectRatio->0.3, MaxRecursion->8, ExclusionsStyle->Directive[Gray,Opacity[0.4]]*)
(*]*)


(* denominator of SpecB must be zero to have a pole \[DoubleLeftRightArrow] roots of eq2d/eq4d *)
eq2d[\[CapitalDelta]ext_,\[Lambda]_] := eq2d[\[CapitalDelta]ext,\[Lambda]] = 1/\[Lambda] + FullSimplify[2 B2d[\[CapitalDelta]ext,\[CapitalDelta]]]
eq4d[\[CapitalDelta]ext_,\[Lambda]_] := eq4d[\[CapitalDelta]ext,\[Lambda]] = 1/\[Lambda] + FullSimplify[2 B4d[\[CapitalDelta]ext,\[CapitalDelta]]]


(* ::Code:: *)
(*eq2d[9/7,10]*)
(*eq4d[9/4,10]*)


Nmax2d = 500;
(* BEWARE: one needs to provide exact values of \[CapitalDelta]ext,\[Lambda], or ones with high enough precision,
	such that FindRoot with WorkingPrecision->100 can work (for example, \[Lambda]=10.81`110) *)
poleSinglet2d[\[CapitalDelta]ext_,\[Lambda]_,n_] := poleSinglet2d[\[CapitalDelta]ext,\[Lambda],n] =
	With[
		{eq = eq2d[\[CapitalDelta]ext,\[Lambda]], \[CapitalDelta]0 = 2\[CapitalDelta]ext, \[Delta] = 10^(-16)},
		(* look for pole in-between 2 MFT poles *)
		\[CapitalDelta] /. FindRoot[eq, {\[CapitalDelta], \[CapitalDelta]0+2n+\[Delta], \[CapitalDelta]0+2n+2-\[Delta]}, Method->"Brent", WorkingPrecision->100]
	]
polesList2d[\[CapitalDelta]ext_,\[Lambda]_] := polesList2d[\[CapitalDelta]ext,\[Lambda]] = ParallelTable[poleSinglet2d[\[CapitalDelta]ext,\[Lambda],n], {n,0,Nmax2d}]


(* ::Code:: *)
(*polesList2d[9/7,100] // N // EchoTiming*)
(*(* show differences from the MFT values, that is the anomalous dimensions *)*)
(*% - Table[2(9/7)+2n,{n,0,Nmax2d}] // N // EchoTiming*)


Nmax4d = 300;
(* TODO: need to take care for extra root appearing for large enough \[Lambda] 
      for now, we simply ignore it *)
poleSinglet4d[\[CapitalDelta]ext_,\[Lambda]_,n_] := poleSinglet4d[\[CapitalDelta]ext,\[Lambda],n] =
	With[
		{eq = eq4d[\[CapitalDelta]ext,\[Lambda]], \[CapitalDelta]0 = 2\[CapitalDelta]ext, \[Delta] = 10^(-16)},
		(*look for pole in-between 2 MFT poles*)
		\[CapitalDelta] /. FindRoot[eq, {\[CapitalDelta], \[CapitalDelta]0+2n+\[Delta], \[CapitalDelta]0+2n+2-\[Delta]}, Method->"Brent", WorkingPrecision->100]
	]
polesList4d[\[CapitalDelta]ext_,\[Lambda]_] := polesList4d[\[CapitalDelta]ext,\[Lambda]] = ParallelTable[poleSinglet4d[\[CapitalDelta]ext,\[Lambda],n], {n,0,Nmax4d}]


(* ::Code:: *)
(*polesList4d[E^(7/10),10];*)
(*%[[-10;;]]*)
(*% - Table[2 E^(7/10)+2n,{n,Nmax4d-9,Nmax4d}]*)


(* ::Subsection:: *)
(*(minus) residues of  the singlet spectral function \[DoubleLeftRightArrow] OPE^2 coefficients of singlet operators*)


(* ::Code:: *)
(*(* for example, in d=2 we have the following *)*)
(*Simplify[(1/D[2 B2d[\[CapitalDelta]ext,\[CapitalDelta]],\[CapitalDelta]]) spec\[CapitalGamma]2d[\[CapitalDelta]ext,\[CapitalDelta]]]*)


opeSqList2d[{nstart_,nend_},\[CapitalDelta]ext_,\[Lambda]_] := opeSqList2d[{nstart,nend},\[CapitalDelta]ext,\[Lambda]] =
	With[
		{opeSqFormula = Simplify[(1/D[2 B2d[\[CapitalDelta]ext,\[CapitalDelta]],\[CapitalDelta]]) Spec\[CapitalGamma]2d[\[CapitalDelta]ext,\[CapitalDelta]]], poles = polesList2d[\[CapitalDelta]ext,\[Lambda]][[nstart;;nend]]},
		Parallelize[(opeSqFormula/.{\[CapitalDelta]->#}) & /@ poles]
	]
opeSqList2d[nmax_,args___] := opeSqList2d[{1,nmax},args]


(* ::Code:: *)
(*opeSqList2d[500,9/7,10] // EchoTiming;*)
(*opeSqList2d[{490,500},9/7,10] // EchoTiming;*)
(*N[%,5]*)
(*%%%[[490;;500]]-%%*)


opeSqList4d[{nstart_,nend_},\[CapitalDelta]ext_,\[Lambda]_] := opeSqList4d[{nstart,nend},\[CapitalDelta]ext,\[Lambda]] =
	With[
		{opeSqFormula = Simplify[(1/D[2 B4d[\[CapitalDelta]ext,\[CapitalDelta]],\[CapitalDelta]]) Spec\[CapitalGamma]4d[\[CapitalDelta]ext,\[CapitalDelta]]], poles = polesList4d[\[CapitalDelta]ext,\[Lambda]][[nstart;;nend]]},
		Parallelize[(opeSqFormula/.{\[CapitalDelta]->#}) & /@ poles]
	]
opeSqList4d[nmax_,args___] := opeSqList4d[{1,nmax},args]


(* ::Code:: *)
(*opeSqList4d[300,9/4,10]// EchoTiming;*)
(*N[%,5]*)


(* ::Section::Closed:: *)
(*Epilog*)


End[]; (* End `Private` Context. *)
(* Protect[Evaluate[Context[] <> "*"]]; (* Protect all public symbols in the package. *) *)
EndPackage[];  (* End package Context. *)
