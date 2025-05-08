(* ::Package:: *)

(* ::Section:: *)
(*Prolog*)


BeginPackage["ONinAdS`CFT`"];


(* expose these variables/functions to public *)
S; K; volS; volSO; volSO1; nPW;
SpecDisc; c2;
\[CapitalOmega]; \[CapitalOmega]holdLeft; \[CapitalOmega]holdRight; \[ScriptCapitalB]2d; \[ScriptCapitalB]2dHold; \[ScriptCapitalB]4d; \[ScriptCapitalB]4dHold;
Ffunction; \[Gamma]tBlock2dHold; \[Gamma]tBlock2d; \[Gamma]tBlock4dHold; \[Gamma]tBlock4d;


Begin["`Private`"];


(* ::Title:: *)
(*CFT generalities*)


(* ::Section:: *)
(*Normalizations appearing in Partial Waves/Conformal Blocks*)


(* (2.16)[d-dimensional SYK, AdS loops, and 6j symbols - Liu, Perlmutter, Rosenhaus, Simmons-Duffin] *)
S[d_,{\[CapitalDelta]1_,\[CapitalDelta]2_},\[CapitalDelta]_,J_] := With[
	{\[CapitalDelta]s=d-\[CapitalDelta], \[CapitalDelta]12=\[CapitalDelta]1-\[CapitalDelta]2},
	((\[Pi]^(d/2) Gamma[\[CapitalDelta]-d/2]Gamma[\[CapitalDelta]+J-1]Gamma[(\[CapitalDelta]s+\[CapitalDelta]12+J)/2]Gamma[(\[CapitalDelta]s-\[CapitalDelta]12+J)/2])/(Gamma[\[CapitalDelta]-1]Gamma[d-\[CapitalDelta]+J]Gamma[(\[CapitalDelta]+\[CapitalDelta]12+J)/2]Gamma[(\[CapitalDelta]-\[CapitalDelta]12+J)/2]))
]
K[d_,{\[CapitalDelta]1_,\[CapitalDelta]2_},\[CapitalDelta]_,J_] := (-(1/2))^J S[d,{\[CapitalDelta]1,\[CapitalDelta]2},\[CapitalDelta],J]
K[d_,\[CapitalDelta]_,J_] := K[d,{0,0},\[CapitalDelta],J]


(* ::Code:: *)
(*K[d,\[CapitalDelta],J] // FullSimplify*)
(*K[2,\[CapitalDelta],J] // FullSimplify*)


(* ::Code:: *)
(*(* compare with slightly different formula in d=2 *)*)
(*(* (2.18)[d-dimensional SYK, AdS loops, and 6j symbols - Liu, Perlmutter, Rosenhaus, Simmons-Duffin] *)*)
(*K2d = With[*)
(*	{h=(\[CapitalDelta]+J)/2,hb=(\[CapitalDelta]-J)/2,sh=1-(\[CapitalDelta]+J)/2,shb=1-(\[CapitalDelta]-J)/2,h1=\[CapitalDelta]1/2,hb1=1-\[CapitalDelta]1/2,h2=\[CapitalDelta]2/2,hb2=1-\[CapitalDelta]2/2},*)
(*	\[Pi] Gamma[1-2sh]/Gamma[2shb] (Gamma[sh-(h1-h2)]Gamma[shb-(hb1-hb2)])/(Gamma[h-(h1-h2)]Gamma[hb-(hb1-hb2)])*)
(*]*)
(*(* d=2 case is different by factor of 2^J compared to general d>2 -> CPWs in d=2 have extra factor of 2^J *)*)
(*FullSimplify[K2d-2^J K[2,{\[CapitalDelta]1,\[CapitalDelta]2},\[CapitalDelta],J], Assumptions->{J\[Element]Integers}]*)


(* ::Code:: *)
(*(* following combination doesn't depend on {\[CapitalDelta]1,\[CapitalDelta]2} *)*)
(*K[d,{\[CapitalDelta]1,\[CapitalDelta]2},d-\[CapitalDelta],J]K[d,{d-\[CapitalDelta]1,d-\[CapitalDelta]2},\[CapitalDelta],J] // Simplify*)


(* volumes *)
volS[n_] := (2 \[Pi]^((n+1)/2))/Gamma[(n+1)/2]
volSO[1] = volSO1;
volSO[n_?IntegerQ] := volS[n-1]volSO[n-1] /; n>0


(* ::Code:: *)
(*Table[volS[n],{n,0,3}]*)
(*Table[volSO[n],{n,1,4}]*)


(* it doesn't depend on \[CapitalDelta]1,\[CapitalDelta]2 [Simmons-Duffin, Stanford, Witten] *)
nPWSD[d_,\[CapitalDelta]_,J_] := (K[d,d-\[CapitalDelta],J]K[d,\[CapitalDelta],J]volS[d-2])/volSO[d-1] ((2J+d-2)\[Pi] Gamma[J+1]Gamma[J+d-2])/(2^(d-2) Gamma[J+d/2]^2)
(* (2.35)[d-dimensional SYK, AdS loops, and 6j symbols - Liu, Perlmutter, Rosenhaus, Simmons-Duffin] *)
(* extra factor 1/2^d compared to [Harmonic Analysis and Mean Field Theory - Karateev, Kravchuk, Simmons-Duffin] *)
nPW[d_,\[CapitalDelta]_,J_] := (K[d,d-\[CapitalDelta],J]K[d,\[CapitalDelta],J]volS[d-2])/(2^d volSO[d-1]) ((2J+d-2)\[Pi] Gamma[J+1]Gamma[J+d-2])/(2^(d-2) Gamma[J+d/2]^2)


(* ::Code:: *)
(*nPW[d,\[CapitalDelta],J] // Simplify*)


(* ::Code:: *)
(*nPW[2,\[CapitalDelta],J] // FullSimplify*)
(*nPW[4,\[CapitalDelta],J] // FullSimplify*)


(* ::Section:: *)
(*Spectral function of t-channel identity (t-channel disconnected contribution in MFT)*)


(* s-channel spectral function of t-channel disconnected correlator (MFT) *)
(* (2.53),(3.118)[Harmonic Analysis and Mean Field Theory - Karateev, Kravchuk, Simmons-Duffin] *)
SpecDisc[d_,{\[CapitalDelta]1_,\[CapitalDelta]2_},\[CapitalDelta]_,J_] := 1/S[d,{\[CapitalDelta]1,\[CapitalDelta]2},d-\[CapitalDelta],J] (
	(2^(J-1) Gamma[\[CapitalDelta]-1]Gamma[d/2-\[CapitalDelta]1]Gamma[d/2-\[CapitalDelta]2]Gamma[d/2+J]Gamma[d+J-\[CapitalDelta]])
	/(Gamma[\[CapitalDelta]1]Gamma[\[CapitalDelta]2]Gamma[J+1]Gamma[\[CapitalDelta]-d/2]Gamma[\[CapitalDelta]+J-1])
	(Gamma[(\[CapitalDelta]+J+\[CapitalDelta]1-\[CapitalDelta]2)/2]Gamma[(\[CapitalDelta]+J-\[CapitalDelta]1+\[CapitalDelta]2)/2]Gamma[(J-\[CapitalDelta]+\[CapitalDelta]1+\[CapitalDelta]2)/2]Gamma[(J+\[CapitalDelta]-d+\[CapitalDelta]1+\[CapitalDelta]2)/2])
	/(Gamma[(d-\[CapitalDelta]+J+\[CapitalDelta]1-\[CapitalDelta]2)/2]Gamma[(d-\[CapitalDelta]+J-\[CapitalDelta]1+\[CapitalDelta]2)/2]Gamma[(2d+J-\[CapitalDelta]-\[CapitalDelta]1-\[CapitalDelta]2)/2]Gamma[(d+J+\[CapitalDelta]-\[CapitalDelta]1-\[CapitalDelta]2)/2])
);
SpecDisc[d_,\[CapitalDelta]ext_,\[CapitalDelta]_,J_] := SpecDisc[d,{\[CapitalDelta]ext,\[CapitalDelta]ext},\[CapitalDelta],J]


(* ::Code:: *)
(*SpecDisc[d,\[CapitalDelta]ext,\[CapitalDelta],J]*)


(* ::Code:: *)
(*(* it is shadow-symmetric, as it should *)*)
(*SpecDisc[d,\[CapitalDelta]ext,\[CapitalDelta],J] - (SpecDisc[d,\[CapitalDelta]ext,\[CapitalDelta],J] /. \[CapitalDelta]->d-\[CapitalDelta]) // FullSimplify*)


(* (4.21) [Carmi, di Pietro, Komatsu] *)
c2[d_,\[CapitalDelta]ext_,n_,J_] := (
	((-1)^J (Pochhammer[\[CapitalDelta]ext-d/2+1,n]Pochhammer[\[CapitalDelta]ext,J+n])^2)
	/(J!n!Pochhammer[J+d/2,n]Pochhammer[2\[CapitalDelta]ext+n-d+1,n]Pochhammer[2\[CapitalDelta]ext+2n+J-1,J]Pochhammer[2\[CapitalDelta]ext+n+J-d/2,n])
);


(* ::Code:: *)
(*c2[d,\[CapitalDelta]ext,n,J]//FunctionExpand*)


(* ::Code:: *)
(*(* compare residue of K*SpecDisc with OPEsquared *)*)
(*$Assumptions={\[CapitalDelta]\[Element]Reals,\[CapitalDelta]ext\[Element]Reals,J\[Element]Integers,J>=0,n\[Element]Integers,n>=0};*)
(*- Residue[K[d,d-\[CapitalDelta],J]SpecDisc[d,\[CapitalDelta]ext,\[CapitalDelta],J],{\[CapitalDelta],2\[CapitalDelta]ext+2n+J}]*)
(*% - c2[d,\[CapitalDelta]ext,n,J] // FullSimplify*)


(* ::Section:: *)
(*Building blocks of the 6-j symbol (crossing kernel)*)


(* central function appearing in 6-j symbols *)
(* (3.38)[d-dimensional SYK, AdS loops, and 6j symbols - Liu, Perlmutter, Rosenhaus, Simmons-Duffin] *)
\[CapitalOmega][h1_,h2_,h3_,h4_,h_,hp_,p_] := (
	(Gamma[2h] Gamma[hp-p+1] Gamma[hp-(h1-h2)+(h3-h4)-p+1] Gamma[-hp+(h1-h2)+h+p-1])
	/(Gamma[(h1-h2)+h] Gamma[(h3-h4)+h] Gamma[hp-(h1-h2)+h-p+1])
	HypergeometricPFQ[{hp+(h2-h3),hp-(h1-h4),hp-(h1-h2)+(h3-h4)-p+1,hp-p+1}, {2hp,hp-(h1-h2)+h-p+1,hp-(h1-h2)-h-p+2}, 1]
	+
	(Gamma[2hp] Gamma[hp-(h1-h2)-h-p+1] Gamma[(h1-h3)+h+p-1] Gamma[(h4-h2)+h+p-1])
	/(Gamma[hp+(h2-h3)] Gamma[hp-(h1-h4)] Gamma[hp+(h1-h2)+h+p-1])
	HypergeometricPFQ[{(h1-h3)+h+p-1,(h4-h2)+h+p-1,(h3-h4)+h,(h1-h2)+h}, {hp+(h1-h2)+h+p-1,2h,-hp+(h1-h2)+h+p}, 1]
)
(* for the case of equal external dimensions, we can subtitute 0s *)
\[CapitalOmega][h_,hp_,p_]:=\[CapitalOmega][0,0,0,0,h,hp,p]
\[CapitalOmega]holdLeft[h1_,h1_,h1_,h1_,h_,hp_,p_]:=\[CapitalOmega]holdLeft[h,hp,p]
\[CapitalOmega]holdRight[h1_,h1_,h1_,h1_,h_,hp_,p_]:=\[CapitalOmega]holdRight[h,hp,p]


(* ::Code:: *)
(*\[CapitalOmega][h,hp,p] // Simplify // TraditionalForm*)


(* ::Code:: *)
(*(* for equal external dimensions, second term is just (h<->hp,p<->2-p) *)*)
(*First@\[CapitalOmega][h,hp,p]-(Last@\[CapitalOmega][h,hp,p]/.{h->hp,hp->h,p->2-p})//Simplify*)


(* ::Subsection:: *)
(*d=2*)


(* "building block" of 6j-symbol, general representation in d=2 with quantum numbers (h,hb) *)
(* (3.36)[d-dimensional SYK, AdS loops, and 6j symbols - Liu, Perlmutter, Rosenhaus, Simmons-Duffin] *)
(* extra (-1)^J=(-1)^(h-hb) since our t-channel/u-channel is swapped compared to (3.5) and below in [d-dimensional SYK, AdS loops, and 6j symbols] *)
\[ScriptCapitalB]2dHold[{h1_,hb1_},{h2_,hb2_},{h3_,hb3_},{h4_,hb4_},{h_,hb_},{hp_,hbp_}] :=
	With[
		{sh=1-h,shb=1-hb,sh3=1-h3,shb3=1-hb3,sh4=1-h4,shb4=1-hb4,
		 h12=h1-h2,hb12=hb1-hb2,h34=h3-h4,hb34=hb3-hb4},
		( 
			((-1)^(h-hb)) (-1)^(h-hb-sh3+shb3-sh4+shb4)/4 *
			(Gamma[h+h12] Gamma[h-h12] Gamma[shb-hb34] Gamma[shb+hb34])/(Gamma[2h]Gamma[2-2hb])*
			Sin[\[Pi](hbp-hb1-hb4)] Sin[\[Pi](hbp-hb2-hb3)]*
			\[CapitalOmega]holdLeft[h1,h2,h3,h4,h,hp,h2+h3] \[CapitalOmega]holdRight[hb1,hb2,hb3,hb4,shb,hbp,hb2+hb3]
		)
	]
\[ScriptCapitalB]2d[{h1_,hb1_},{h2_,hb2_},{h3_,hb3_},{h4_,hb4_},{h_,hb_},{hp_,hbp_}] :=
	\[ScriptCapitalB]2dHold[{h1,hb1},{h2,hb2},{h3,hb3},{h4,hb4},{h,hb},{hp,hbp}] /. {\[CapitalOmega]holdLeft->\[CapitalOmega],\[CapitalOmega]holdRight->\[CapitalOmega]}


(* external scalars with dimensions \[CapitalDelta]ext, and exchanged symmetric traceless tensor (STT) with (\[CapitalDelta]p,Jp),
   which for nonzero Jp consists of (hp,hbp)+(hbp,hp) *)
(* version without substituting \[CapitalOmega] *)
\[ScriptCapitalB]2dHold[\[CapitalDelta]ext_,{\[CapitalDelta]_,J_},{\[CapitalDelta]p_,Jp_}] :=
	With[
		{hext={\[CapitalDelta]ext/2,\[CapitalDelta]ext/2}, h=(\[CapitalDelta]+J)/2, hb=(\[CapitalDelta]-J)/2, hp=(\[CapitalDelta]p+Jp)/2, hbp=(\[CapitalDelta]p-Jp)/2},
		(1/(1+KroneckerDelta[0,Jp]))*
		(\[ScriptCapitalB]2dHold[hext,hext,hext,hext,{h,hb},{hp,hbp}]+\[ScriptCapitalB]2dHold[hext,hext,hext,hext,{h,hb},{hbp,hp}])
	]
(* version with substituting \[CapitalOmega] *)
\[ScriptCapitalB]2d[\[CapitalDelta]ext_,{\[CapitalDelta]_,J_},{\[CapitalDelta]p_,Jp_}] := \[ScriptCapitalB]2dHold[\[CapitalDelta]ext,{\[CapitalDelta],J},{\[CapitalDelta]p,Jp}] /. {\[CapitalOmega]holdLeft->\[CapitalOmega],\[CapitalOmega]holdRight->\[CapitalOmega]}


(* ::Subsection:: *)
(*d=4*)


(* (2.38) and (2.39) [d-dimensional SYK, AdS loops, and 6j symbols - Liu, Perlmutter, Rosenhaus, Simmons-Duffin] *)
t0[d_,J_] := 1/(2^d volSO[d-1]) ((Gamma[(d-2)/2] Gamma[J+d-2])/(2^J Gamma[d-2] Gamma[J+(d-2)/2]))
(* (3.14) [d-dimensional SYK, AdS loops, and 6j symbols - Liu, Perlmutter, Rosenhaus, Simmons-Duffin] *)
\[Alpha][d_,{\[CapitalDelta]1_,\[CapitalDelta]2_,\[CapitalDelta]3_,\[CapitalDelta]4_},\[CapitalDelta]_,J_] := With[
	{s\[CapitalDelta]=d-\[CapitalDelta],\[CapitalDelta]12=\[CapitalDelta]1-\[CapitalDelta]2,\[CapitalDelta]34=\[CapitalDelta]3-\[CapitalDelta]4},
	-(t0[d,J]/2^(d+1)) (2\[Pi])^(d-2) (Gamma[J+1]Gamma[\[CapitalDelta]-d/2])/(Gamma[J+d/2]Gamma[\[CapitalDelta]-1]) (Gamma[(\[CapitalDelta]12+J+\[CapitalDelta])/2]Gamma[(-\[CapitalDelta]12+J+\[CapitalDelta])/2]Gamma[(\[CapitalDelta]34+J+s\[CapitalDelta])/2]Gamma[(-\[CapitalDelta]34+J+s\[CapitalDelta])/2])/(Gamma[J+\[CapitalDelta]]Gamma[J+d-\[CapitalDelta]])
 ]
(* (3.43) [d-dimensional SYK, AdS loops, and 6j symbols - Liu, Perlmutter, Rosenhaus, Simmons-Duffin] *)
\[CapitalTheta][{\[CapitalDelta]1_,\[CapitalDelta]2_,\[CapitalDelta]3_,\[CapitalDelta]4_},x_] := (4\[Pi]^2)/(Gamma[(\[CapitalDelta]3+\[CapitalDelta]2-x)/2]Gamma[1-(\[CapitalDelta]3+\[CapitalDelta]2-x)/2]Gamma[(\[CapitalDelta]1+\[CapitalDelta]4-x)/2]Gamma[1-(\[CapitalDelta]1+\[CapitalDelta]4-x)/2])


(* "building block" of 6j-symbol, d=4 *)
(* version without substituting \[CapitalOmega] (3.42) [d-dimensional SYK, AdS loops, and 6j symbols - Liu, Perlmutter, Rosenhaus, Simmons-Duffin] *)
(* extra (-1)^J since our t-channel/u-channel is swapped compared to (3.5) and below in [d-dimensional SYK, AdS loops, and 6j symbols] *)
(* CARE: make sure to make \[CapitalOmega]holdRight always the one with s\[CapitalDelta], such that we pick all residues (even with \[CapitalOmega]holdLeft not expanded) *)
(* version without substituting \[CapitalOmega] *)
\[ScriptCapitalB]4dHold[{\[CapitalDelta]1_,\[CapitalDelta]2_,\[CapitalDelta]3_,\[CapitalDelta]4_},{\[CapitalDelta]_,J_},{\[CapitalDelta]p_,Jp_}] :=
	With[
		{s\[CapitalDelta]=4-\[CapitalDelta],sJp=-2-Jp},
		(-1)^J (-1)^(J) \[Alpha][4,{\[CapitalDelta]1,\[CapitalDelta]2,\[CapitalDelta]3,\[CapitalDelta]4},\[CapitalDelta],J]( 
			 \[CapitalTheta][{\[CapitalDelta]1,\[CapitalDelta]2,\[CapitalDelta]3,\[CapitalDelta]4},\[CapitalDelta]p+sJp]  \[CapitalOmega]holdLeft[\[CapitalDelta]1/2,\[CapitalDelta]2/2,\[CapitalDelta]3/2,\[CapitalDelta]4/2,( \[CapitalDelta]+J)/2,(\[CapitalDelta]p+ Jp)/2,(\[CapitalDelta]2+\[CapitalDelta]3)/2-1] \[CapitalOmega]holdRight[\[CapitalDelta]1/2,\[CapitalDelta]2/2,\[CapitalDelta]3/2,\[CapitalDelta]4/2,(s\[CapitalDelta]+J)/2,(\[CapitalDelta]p+sJp)/2,(\[CapitalDelta]2+\[CapitalDelta]3)/2-1]
			-\[CapitalTheta][{\[CapitalDelta]1,\[CapitalDelta]2,\[CapitalDelta]3,\[CapitalDelta]4},\[CapitalDelta]p+ Jp]  \[CapitalOmega]holdLeft[\[CapitalDelta]1/2,\[CapitalDelta]2/2,\[CapitalDelta]3/2,\[CapitalDelta]4/2,( \[CapitalDelta]+J)/2,(\[CapitalDelta]p+sJp)/2,(\[CapitalDelta]2+\[CapitalDelta]3)/2-1] \[CapitalOmega]holdRight[\[CapitalDelta]1/2,\[CapitalDelta]2/2,\[CapitalDelta]3/2,\[CapitalDelta]4/2,(s\[CapitalDelta]+J)/2,(\[CapitalDelta]p+ Jp)/2,(\[CapitalDelta]2+\[CapitalDelta]3)/2-1] (*Jp<->sJp*)
			-\[CapitalTheta][{\[CapitalDelta]1,\[CapitalDelta]2,\[CapitalDelta]3,\[CapitalDelta]4},\[CapitalDelta]p+sJp] \[CapitalOmega]holdRight[\[CapitalDelta]1/2,\[CapitalDelta]2/2,\[CapitalDelta]3/2,\[CapitalDelta]4/2,(s\[CapitalDelta]+J)/2,(\[CapitalDelta]p+ Jp)/2,(\[CapitalDelta]2+\[CapitalDelta]3)/2-1]  \[CapitalOmega]holdLeft[\[CapitalDelta]1/2,\[CapitalDelta]2/2,\[CapitalDelta]3/2,\[CapitalDelta]4/2,( \[CapitalDelta]+J)/2,(\[CapitalDelta]p+sJp)/2,(\[CapitalDelta]2+\[CapitalDelta]3)/2-1] (*\[CapitalDelta]<->s\[CapitalDelta]*)
			+\[CapitalTheta][{\[CapitalDelta]1,\[CapitalDelta]2,\[CapitalDelta]3,\[CapitalDelta]4},\[CapitalDelta]p+ Jp] \[CapitalOmega]holdRight[\[CapitalDelta]1/2,\[CapitalDelta]2/2,\[CapitalDelta]3/2,\[CapitalDelta]4/2,(s\[CapitalDelta]+J)/2,(\[CapitalDelta]p+sJp)/2,(\[CapitalDelta]2+\[CapitalDelta]3)/2-1]  \[CapitalOmega]holdLeft[\[CapitalDelta]1/2,\[CapitalDelta]2/2,\[CapitalDelta]3/2,\[CapitalDelta]4/2,( \[CapitalDelta]+J)/2,(\[CapitalDelta]p+ Jp)/2,(\[CapitalDelta]2+\[CapitalDelta]3)/2-1] (*Jp<->sJp,\[CapitalDelta]<->s\[CapitalDelta]*)
		)
	]
(* version with substituting \[CapitalOmega] *)
\[ScriptCapitalB]4d[{\[CapitalDelta]1_,\[CapitalDelta]2_,\[CapitalDelta]3_,\[CapitalDelta]4_},{\[CapitalDelta]_,J_},{\[CapitalDelta]p_,Jp_}] := \[ScriptCapitalB]4dHold[{\[CapitalDelta]1,\[CapitalDelta]2,\[CapitalDelta]3,\[CapitalDelta]4},{\[CapitalDelta],J},{\[CapitalDelta]p,Jp}] /. {\[CapitalOmega]holdLeft->\[CapitalOmega],\[CapitalOmega]holdRight->\[CapitalOmega]}
(* version with scalars on external legs with equal scaling dimension *)
\[ScriptCapitalB]4dHold[\[CapitalDelta]ext_,{\[CapitalDelta]_,J_},{\[CapitalDelta]p_,Jp_}] := \[ScriptCapitalB]4dHold[{\[CapitalDelta]ext,\[CapitalDelta]ext,\[CapitalDelta]ext,\[CapitalDelta]ext},{\[CapitalDelta],J},{\[CapitalDelta]p,Jp}]
\[ScriptCapitalB]4d[\[CapitalDelta]ext_,{\[CapitalDelta]_,J_},{\[CapitalDelta]p_,Jp_}] := \[ScriptCapitalB]4dHold[\[CapitalDelta]ext,{\[CapitalDelta],J},{\[CapitalDelta]p,Jp}] /. {\[CapitalOmega]holdLeft->\[CapitalOmega],\[CapitalOmega]holdRight->\[CapitalOmega]}


(* ::Section:: *)
(*Calculation  \[LongDash] contribution of a single t-channel conformal block to anomalous dimensions*)


(* ::Subsection:: *)
(*d=2*)


(* ::Code:: *)
(*$Assumptions = {\[CapitalDelta]\[Element]Reals,\[CapitalDelta]ext\[Element]Reals,J\[Element]Integers,Jp\[Element]Integers,J>=0,Jp>=0,n\[Element]Integers,n>=0};*)
(*(* need to be careful with substituting {\[CapitalDelta]p,Jp}->{0,0} *)*)
(*\[ScriptCapitalB]2d00 = (\[ScriptCapitalB]2d[\[CapitalDelta]ext,{\[CapitalDelta],J},{\[CapitalDelta]p,0}]//Simplify//FunctionExpand)/.\[CapitalDelta]p->0 // Simplify*)
(*SpecDisc2d = SpecDisc[2,{\[CapitalDelta]ext,\[CapitalDelta]ext},\[CapitalDelta],J] // Simplify;*)
(*(2^(1-J) \[ScriptCapitalB]2d00/nPW[2,\[CapitalDelta],J])/(SpecDisc2d) // FullSimplify*)
(*(* the t-channel identity spectral function is same as \[ScriptCapitalB]2d00 (up to normalization) *)*)
(*(* thus, the calculation via CrossKernel/SpecDisc agrees with B/B00 *)*)
(*(* Additional comments:*)
(*   - factor volSO1 appears most likely because it was set to 1 in (2.38);*)
(*   - factor of 2^(1-J) counteracting a different normalization of CPW in d=2*)
(*         (see Section 2 of [Simmons-Duffin, Stanford, Witten] for more details) *)*)


(* ::Code:: *)
(*(* contribution to the anomalous dimensions is given by *)*)
(*CrossKernel2d[\[CapitalDelta]ext_,{\[CapitalDelta]_,J_},{\[CapitalDelta]p_,Jp_}] := (2^(1-J) \[ScriptCapitalB]2dHold[\[CapitalDelta]ext,{\[CapitalDelta],J},{\[CapitalDelta]p,Jp}])/nPW[2,\[CapitalDelta],J] /. {\[CapitalOmega]holdRight->\[CapitalOmega],volSO1->1}*)
(*\[Gamma]tContrib2d = Residue[CrossKernel2d[\[CapitalDelta]ext,{\[CapitalDelta],J},{\[CapitalDelta]p,Jp}]/SpecDisc2d,{\[CapitalDelta],2\[CapitalDelta]ext+2n+J}] // Simplify*)


(* ::Code:: *)
(*(* check, that this result is indeed the one defined in the next subsection *)*)
(*\[Gamma]tContrib2d/\[Gamma]tBlock2dHold[\[CapitalDelta]ext,{n,J},{\[CapitalDelta]p,Jp}] // FullSimplify*)


(* ::Subsection:: *)
(*d=4*)


(* ::Code:: *)
(*$Assumptions = {\[CapitalDelta]\[Element]Reals,\[CapitalDelta]ext\[Element]Reals,J\[Element]Integers,Jp\[Element]Integers,J>=0,Jp>=0,n\[Element]Integers,n>=0};*)
(*(* need to be careful with substituting {\[CapitalDelta]p,Jp}->{0,0} *)*)
(*\[ScriptCapitalB]4d00 = (\[ScriptCapitalB]4d[\[CapitalDelta]ext,{\[CapitalDelta],J},{\[CapitalDelta]p,0}]//Simplify//FunctionExpand)/.\[CapitalDelta]p->0 // FullSimplify*)
(*SpecDisc4d = SpecDisc[4,{\[CapitalDelta]ext,\[CapitalDelta]ext},\[CapitalDelta],J] // Simplify;*)
(*(* the t-channel identity spectral function is same as \[ScriptCapitalB]4d00 (up to normalization) *)*)
(*(* thus, the calculation via CrossKernel/SpecDisc agrees with B/B00 *)*)
(*(\[ScriptCapitalB]4d00/nPW[4,\[CapitalDelta],J])/(SpecDisc4d) // FullSimplify*)


(* ::Code:: *)
(*(* contribution to the anomalous dimensions is given by *)*)
(*CrossKernel4d[\[CapitalDelta]ext_,{\[CapitalDelta]_,J_},{\[CapitalDelta]p_,Jp_}] := \[ScriptCapitalB]4dHold[\[CapitalDelta]ext,{\[CapitalDelta],J},{\[CapitalDelta]p,Jp}]/nPW[4,\[CapitalDelta],J] /. \[CapitalOmega]holdRight->\[CapitalOmega]*)
(*\[Gamma]tContrib4d = Residue[CrossKernel4d[\[CapitalDelta]ext,{\[CapitalDelta],J},{\[CapitalDelta]p,Jp}]/SpecDisc4d,{\[CapitalDelta],2\[CapitalDelta]ext+2n+J}] // Simplify*)


(* ::Code:: *)
(*(* check, that this result is indeed the one defined in the next subsection *)*)
(*\[Gamma]tContrib4d/\[Gamma]tBlock4dHold[\[CapitalDelta]ext,{n,J},{\[CapitalDelta]p,Jp}] // FullSimplify*)
(*(* need to help Mathematica with the final comparison *)*)
(*Table[%, {Jp,0,10}] // FullSimplify*)


(* ::Subsection:: *)
(*Compare Asymptotics for leading twist*)


(* ::Code:: *)
(*(* quick check - the leading-twist large J asymptotics agree in d=2 and d=4 *)*)
(*ratio2dand4dLeadingTwist=\[Gamma]tContrib4d/\[Gamma]tContrib2d /.{n->0} //FullSimplify*)
(*N[ParallelTable[{10^logJ,ratio2dand4dLeadingTwist/.\[CapitalOmega]holdLeft->\[CapitalOmega]}/.{\[CapitalDelta]ext->E,\[CapitalDelta]p->2E+1+0.12`40,J->10^logJ},{Jp,0,2},{logJ,0,3,1}],7]//EchoTiming*)
(*ListPlot[%,ScalingFunctions->{"Log","SignedLog"}]*)


(* ::Section:: *)
(*Result \[LongDash] contribution of a single t-channel conformal block to anomalous dimensions*)


Ffunction[d_,\[CapitalDelta]ext_,n_,a_] := (
	((-1)^n/n!)(Gamma[d-2\[CapitalDelta]ext-n]/Gamma[d-2\[CapitalDelta]ext-2n])(Gamma[d/2-\[CapitalDelta]ext-n]^2/Gamma[d/2-\[CapitalDelta]ext]^2) (Gamma[a+n]/Gamma[a-n])*
	HypergeometricPFQ[{-n,-n,d/2-\[CapitalDelta]ext-n,d/2-\[CapitalDelta]ext-n}, {d-2\[CapitalDelta]ext-2n, 1-a-n, a-n}, 1]
)


(* ::Subsection:: *)
(*d=2*)


(* general formula for contribution of a t-channel conformal block to anomalous dimension of s-channel double-trace operators *)
(* arbitary-twist version of (3.55)[d-dimensional SYK, AdS loops, and 6j symbols - Liu, Perlmutter, Rosenhaus, Simmons-Duffin],
	if one picks just the first \[CapitalOmega]-term *)
\[Gamma]tBlock2dHold[\[CapitalDelta]ext_, {n_,J_}, {\[CapitalDelta]p_, Jp_}] := With[
	{sJp=-Jp},
	-(2 Gamma[\[CapitalDelta]ext]^2 (Gamma[\[CapitalDelta]ext+J+n]^2/Gamma[2(\[CapitalDelta]ext+J+n)])) (Sin[\[Pi](\[CapitalDelta]ext-(\[CapitalDelta]p-Jp)/2)]^2/\[Pi]^2) (Gamma[J+n+1]/Gamma[2\[CapitalDelta]ext+J+n-1])*
	(1/(1+KroneckerDelta[0,Jp]))*(
	  (Gamma[\[CapitalDelta]p+sJp]/Gamma[(\[CapitalDelta]p+sJp)/2]^2) Ffunction[2,\[CapitalDelta]ext,n,(\[CapitalDelta]p+sJp)/2]	\[CapitalOmega]holdLeft[\[CapitalDelta]ext+J+n, (\[CapitalDelta]p+ Jp)/2, \[CapitalDelta]ext]
	 +(Gamma[\[CapitalDelta]p+ Jp]/Gamma[(\[CapitalDelta]p+ Jp)/2]^2) Ffunction[2,\[CapitalDelta]ext,n,(\[CapitalDelta]p+ Jp)/2]	\[CapitalOmega]holdLeft[\[CapitalDelta]ext+J+n, (\[CapitalDelta]p+sJp)/2, \[CapitalDelta]ext]
	)
]
\[Gamma]tBlock2d[\[CapitalDelta]ext_, {n_,J_}, {\[CapitalDelta]p_, Jp_}] := \[Gamma]tBlock2dHold[\[CapitalDelta]ext, {n,J}, {\[CapitalDelta]p, Jp}] /. {\[CapitalOmega]holdLeft->\[CapitalOmega]}
\[Gamma]tBlock2d[\[CapitalDelta]ext_, {n_,J_}, \[CapitalDelta]p_] := \[Gamma]tBlock2d[\[CapitalDelta]ext, {n,J}, {\[CapitalDelta]p, 0}]


(* ::Subsection:: *)
(*d=4*)


(* result after performing substitution in 
	(3.56)[d-dimensional SYK, AdS loops, and 6j symbols - Liu, Perlmutter, Rosenhaus, Simmons-Duffin] *)
\[Gamma]0tBlock4dHold[\[CapitalDelta]ext_,J_,{\[CapitalDelta]p_,Jp_}] := (
	((-1+Cos[\[Pi] (Jp+2 \[CapitalDelta]ext-\[CapitalDelta]p)]) J! Gamma[\[CapitalDelta]ext]^2 Gamma[J+\[CapitalDelta]ext]^2 (
		Gamma[1/2 (-2-Jp+\[CapitalDelta]p)]^2 Gamma[Jp+\[CapitalDelta]p] \[CapitalOmega]holdLeft[\[CapitalDelta]ext/2,\[CapitalDelta]ext/2,\[CapitalDelta]ext/2,\[CapitalDelta]ext/2,J+\[CapitalDelta]ext,1/2 (-2-Jp+\[CapitalDelta]p),-1+\[CapitalDelta]ext]
		-Gamma[-2-Jp+\[CapitalDelta]p] Gamma[(Jp+\[CapitalDelta]p)/2]^2 \[CapitalOmega]holdLeft[\[CapitalDelta]ext/2,\[CapitalDelta]ext/2,\[CapitalDelta]ext/2,\[CapitalDelta]ext/2,J+\[CapitalDelta]ext,(Jp+\[CapitalDelta]p)/2,-1+\[CapitalDelta]ext])
	)/(\[Pi]^2 Gamma[2 (J+\[CapitalDelta]ext)] Gamma[-1+J+2 \[CapitalDelta]ext] Gamma[1/2 (-2-Jp+\[CapitalDelta]p)]^2 Gamma[(Jp+\[CapitalDelta]p)/2]^2)
)
\[Gamma]0tBlock4d[\[CapitalDelta]ext_, J_, {\[CapitalDelta]p_, Jp_}] := \[Gamma]0tBlock4dHold[\[CapitalDelta]ext, J, {\[CapitalDelta]p, Jp}] /. {\[CapitalOmega]holdLeft->\[CapitalOmega]}


(* general formula for contribution of a t-channel conformal block to anomalous dimension of s-channel double-trace operators *)
(* arbitary-twist version of (3.56)[d-dimensional SYK, AdS loops, and 6j symbols - Liu, Perlmutter, Rosenhaus, Simmons-Duffin] *)
\[Gamma]tBlock4dHold[\[CapitalDelta]ext_, {n_,J_}, {\[CapitalDelta]p_, Jp_}] := With[
	{sJp=-2-Jp},
	(2 Gamma[\[CapitalDelta]ext]^2 (Gamma[\[CapitalDelta]ext+J+n]^2/Gamma[2(\[CapitalDelta]ext+J+n)])) * (Sin[\[Pi](\[CapitalDelta]ext-(\[CapitalDelta]p-Jp)/2)]^2/\[Pi]^2) *
	(Gamma[J+n+1]/Gamma[2\[CapitalDelta]ext+J+n-1]) *
	((J+n+1)/(J+1))*((2\[CapitalDelta]ext+J+n-2)/(2\[CapitalDelta]ext+J+2n-2))* ( 
	   (Gamma[\[CapitalDelta]p+sJp]/Gamma[(\[CapitalDelta]p+sJp)/2]^2) Ffunction[4,\[CapitalDelta]ext,n,(\[CapitalDelta]p+sJp)/2]	\[CapitalOmega]holdLeft[\[CapitalDelta]ext+J+n, (\[CapitalDelta]p+Jp)/2, \[CapitalDelta]ext-1]
	  -(Gamma[\[CapitalDelta]p+Jp]/Gamma[(\[CapitalDelta]p+Jp)/2]^2) Ffunction[4,\[CapitalDelta]ext,n,(\[CapitalDelta]p+Jp)/2]	\[CapitalOmega]holdLeft[\[CapitalDelta]ext+J+n, (\[CapitalDelta]p+sJp)/2, \[CapitalDelta]ext-1]
	)
]
\[Gamma]tBlock4d[\[CapitalDelta]ext_, {n_,J_}, {\[CapitalDelta]p_, Jp_}] := \[Gamma]tBlock4dHold[\[CapitalDelta]ext, {n,J}, {\[CapitalDelta]p, Jp}] /. {\[CapitalOmega]holdLeft->\[CapitalOmega]}
\[Gamma]tBlock4d[\[CapitalDelta]ext_, {n_,J_}, \[CapitalDelta]p_] := \[Gamma]tBlock4d[\[CapitalDelta]ext, {n,J}, {\[CapitalDelta]p, 0}]


(* ::Code:: *)
(*(* comparison with the original n=0 formula *)*)
(*\[Gamma]tBlock4dHold[\[CapitalDelta]ext,{0,J},{\[CapitalDelta]p,Jp}]/\[Gamma]0tBlock4dHold[\[CapitalDelta]ext,J,{\[CapitalDelta]p,Jp}] // FullSimplify*)


(* ::Subsection:: *)
(*Check the large J asymptotics*)


(* ::Code:: *)
(*$Assumptions = {\[CapitalDelta]\[Element]Reals,\[CapitalDelta]p\[Element]Reals,\[CapitalDelta]ext\[Element]Reals,\[CapitalDelta]>0,\[CapitalDelta]p>0,\[CapitalDelta]ext>0,J\[Element]Integers,Jp\[Element]Integers,J>=0,Jp>=0,n\[Element]Integers,n>=0,J>Jp/2-2\[CapitalDelta]ext+\[CapitalDelta]p/2};*)
(*asymptoticshouldbe=(2 Gamma[\[CapitalDelta]p+Jp]/(Gamma[(\[CapitalDelta]p+Jp)/2]^2) Gamma[\[CapitalDelta]ext]^2/Gamma[\[CapitalDelta]ext-(\[CapitalDelta]p-Jp)/2]^2 J^(-(\[CapitalDelta]p-Jp)))*)
(*asymptoticratio2d=\[Gamma]tBlock2d[\[CapitalDelta]ext,{0,J},{\[CapitalDelta]p,Jp}]/asymptoticshouldbe //Simplify // FunctionExpand;*)


(* ::Code:: *)
(*table=N[ParallelTable[{10^logJ,asymptoticratio2d}/.{\[CapitalDelta]ext->E,\[CapitalDelta]p->2E+1+0.12`200,J->10^logJ},{Jp,0,3},{logJ,2,3.5,1/2}],7]//EchoTiming*)
(*ListPlot[table,ScalingFunctions->{"Log","SignedLog"}]*)


(* ::Code:: *)
(*$Assumptions = {\[CapitalDelta]\[Element]Reals,\[CapitalDelta]p\[Element]Reals,\[CapitalDelta]ext\[Element]Reals,\[CapitalDelta]>0,\[CapitalDelta]p>0,\[CapitalDelta]ext>0,J\[Element]Integers,Jp\[Element]Integers,J>=0,Jp>=0,n\[Element]Integers,n>=0,J>Jp/2-2\[CapitalDelta]ext+\[CapitalDelta]p/2};*)
(*asymptoticshouldbe=(2 Gamma[\[CapitalDelta]p+Jp]/Gamma[(\[CapitalDelta]p+Jp)/2]^2 Gamma[\[CapitalDelta]ext]^2/Gamma[\[CapitalDelta]ext-(\[CapitalDelta]p-Jp)/2]^2 J^(-(\[CapitalDelta]p-Jp)))*)
(*asymptoticratio4d=\[Gamma]tBlock4d[\[CapitalDelta]ext,{0,J},{\[CapitalDelta]p,Jp}]/asymptoticshouldbe;*)


(* ::Code:: *)
(*table4d=N[ParallelTable[{10^logJ,asymptoticratio4d}/.{\[CapitalDelta]ext->E^(7/4),\[CapitalDelta]p->8.6`200,J->10^logJ},{Jp,0,4},{logJ,2,3.5,1/2}],7]//EchoTiming*)
(*ListPlot[table4d,ScalingFunctions->{"Log","SignedLog"}]*)


(* ::Subsubsection:: *)
(*The "second" \[CapitalOmega]-terms are dominant (those with \[CapitalDelta]p-Jp in the second argument)*)


(* ::Code:: *)
(*\[Gamma]tBlock4dHoldFirst[\[CapitalDelta]ext_, {n_,J_}, {\[CapitalDelta]p_, Jp_}] := With[*)
(*	{sJp=-2-Jp},*)
(*	(2 Gamma[\[CapitalDelta]ext]^2 (Gamma[\[CapitalDelta]ext+J+n]^2/Gamma[2(\[CapitalDelta]ext+J+n)])) * (Sin[\[Pi](\[CapitalDelta]ext-(\[CapitalDelta]p-Jp)/2)]^2/\[Pi]^2) **)
(*	(Gamma[J+n+1]/Gamma[2\[CapitalDelta]ext+J+n-1]) **)
(*	((J+n+1)/(J+1))*((2\[CapitalDelta]ext+J+n-2)/(2\[CapitalDelta]ext+J+2n-2))* ( *)
(*	   (Gamma[\[CapitalDelta]p+sJp]/Gamma[(\[CapitalDelta]p+sJp)/2]^2) Ffunction[4,\[CapitalDelta]ext,n,(\[CapitalDelta]p+sJp)/2]	\[CapitalOmega]holdLeft[\[CapitalDelta]ext+J+n, (\[CapitalDelta]p+Jp)/2, \[CapitalDelta]ext-1]*)
(*	  (*-(Gamma[\[CapitalDelta]p+Jp]/Gamma[(\[CapitalDelta]p+Jp)/2]^2) Ffunction[4,\[CapitalDelta]ext,n,(\[CapitalDelta]p+Jp)/2]	\[CapitalOmega]holdLeft[\[CapitalDelta]ext+J+n, (\[CapitalDelta]p+sJp)/2, \[CapitalDelta]ext-1]*)*)
(*	)*)
(*]*)
(*\[Gamma]tBlock4dFirst[\[CapitalDelta]ext_, {n_,J_}, {\[CapitalDelta]p_, Jp_}] := \[Gamma]tBlock4dHoldFirst[\[CapitalDelta]ext, {n,J}, {\[CapitalDelta]p, Jp}] /. {\[CapitalOmega]holdLeft->\[CapitalOmega]}*)
(*\[Gamma]tBlock4dHoldSecond[\[CapitalDelta]ext_, {n_,J_}, {\[CapitalDelta]p_, Jp_}] := With[*)
(*	{sJp=-2-Jp},*)
(*	(2 Gamma[\[CapitalDelta]ext]^2 (Gamma[\[CapitalDelta]ext+J+n]^2/Gamma[2(\[CapitalDelta]ext+J+n)])) * (Sin[\[Pi](\[CapitalDelta]ext-(\[CapitalDelta]p-Jp)/2)]^2/\[Pi]^2) **)
(*	(Gamma[J+n+1]/Gamma[2\[CapitalDelta]ext+J+n-1]) **)
(*	((J+n+1)/(J+1))*((2\[CapitalDelta]ext+J+n-2)/(2\[CapitalDelta]ext+J+2n-2))* ( *)
(*	   (*(Gamma[\[CapitalDelta]p+sJp]/Gamma[(\[CapitalDelta]p+sJp)/2]^2) Ffunction[4,\[CapitalDelta]ext,n,(\[CapitalDelta]p+sJp)/2]	\[CapitalOmega]holdLeft[\[CapitalDelta]ext+J+n, (\[CapitalDelta]p+Jp)/2, \[CapitalDelta]ext-1]*)*)
(*	  -(Gamma[\[CapitalDelta]p+Jp]/Gamma[(\[CapitalDelta]p+Jp)/2]^2) Ffunction[4,\[CapitalDelta]ext,n,(\[CapitalDelta]p+Jp)/2]	\[CapitalOmega]holdLeft[\[CapitalDelta]ext+J+n, (\[CapitalDelta]p+sJp)/2, \[CapitalDelta]ext-1]*)
(*	)*)
(*]*)
(*\[Gamma]tBlock4dSecond[\[CapitalDelta]ext_, {n_,J_}, {\[CapitalDelta]p_, Jp_}] := \[Gamma]tBlock4dHoldSecond[\[CapitalDelta]ext, {n,J}, {\[CapitalDelta]p, Jp}] /. {\[CapitalOmega]holdLeft->\[CapitalOmega]}*)


(* ::Code:: *)
(*$Assumptions = {\[CapitalDelta]\[Element]Reals,\[CapitalDelta]p\[Element]Reals,\[CapitalDelta]ext\[Element]Reals,\[CapitalDelta]>0,\[CapitalDelta]p>0,\[CapitalDelta]ext>0,J\[Element]Integers,Jp\[Element]Integers,J>=0,Jp>=0,n\[Element]Integers,n>=0,J>Jp/2-2\[CapitalDelta]ext+\[CapitalDelta]p/2};*)
(*asymptoticshouldbe=(2 Gamma[\[CapitalDelta]p+Jp]/Gamma[(\[CapitalDelta]p+Jp)/2]^2 Gamma[\[CapitalDelta]ext]^2/Gamma[\[CapitalDelta]ext-(\[CapitalDelta]p-Jp)/2]^2 J^(-(\[CapitalDelta]p-Jp)))*)
(*asymptoticratioFirst4d=\[Gamma]tBlock4dFirst[\[CapitalDelta]ext,{0,J},{\[CapitalDelta]p,Jp}]/asymptoticshouldbe;*)
(*asymptoticratioSecond4d=\[Gamma]tBlock4dSecond[\[CapitalDelta]ext,{0,J},{\[CapitalDelta]p,Jp}]/asymptoticshouldbe;*)


(* ::Code:: *)
(*JpMax=1; logJMax=3.5;*)
(*table4dFirst=N[ParallelTable[{10^logJ,asymptoticratioFirst4d}/.{\[CapitalDelta]ext->E^(7/4),\[CapitalDelta]p->8.2`200,J->10^logJ},{Jp,0,JpMax},{logJ,2,logJMax,1/2}],7]//EchoTiming*)
(*table4dSecond=N[ParallelTable[{10^logJ,asymptoticratioSecond4d}/.{\[CapitalDelta]ext->E^(7/4),\[CapitalDelta]p->8.2`200,J->10^logJ},{Jp,0,JpMax},{logJ,2,logJMax,1/2}],7]//EchoTiming*)
(*{ListPlot[table4dFirst,ScalingFunctions->{"Log","SignedLog"}],ListPlot[table4dSecond,ScalingFunctions->{"Log","SignedLog"}]}*)


(* ::Section::Closed:: *)
(*Epilog*)


End[]; (* End `Private` Context. *)
(* Protect[Evaluate[Context[] <> "*"]]; (* Protect all public symbols in the package. *) *)
EndPackage[];  (* End package Context. *)
