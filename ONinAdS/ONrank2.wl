(* ::Package:: *)

(* ::Section:: *)
(*Prolog*)


BeginPackage["ONinAdS`ONrank2`", "ONinAdS`CFT`", "ONinAdS`ONsinglet`"];


(* clear all declared variables/functions which memoize their values *)
ClearAll[\[Gamma]tBlockList2d,\[Gamma]tBlockList4d];
ClearAll[\[Gamma]List2d, \[Gamma]List4d];
ClearAll[\[Gamma]2d, \[Gamma]4d];


(* expose these variables/functions to public *)
\[Gamma]tBlockList2d; \[Gamma]tBlockList4d;
\[Gamma]List2d::usage = "\[Gamma]List2d[{tnstart,tnend},\[CapitalDelta]ext,\[Lambda],{n,J}] computes the contributions of t-channel non-MFT scalar operators to the anomalous dimensions of O(N) rank-2 {n,J} operators in d=2";
\[Gamma]List4d::usage = "\[Gamma]List4d[{tnstart,tnend},\[CapitalDelta]ext,\[Lambda],{n,J}] computes the contributions of t-channel non-MFT scalar operators to the anomalous dimensions of O(N) rank-2 {n,J} operators in d=4";
\[Gamma]2d::usage = "\[Gamma]2d[{tnstart,tnend},\[CapitalDelta]ext,\[Lambda],{n,J}] computes the anomalous dimensions of O(N) rank-2 {n,J} double-twist operators in d=2";
\[Gamma]4d::usage = "\[Gamma]4d[{tnstart,tnend},\[CapitalDelta]ext,\[Lambda],{n,J}] computes the anomalous dimensions of O(N) rank-2 {n,J} double-twist operators in d=4";


Begin["`Private`"];


(* ::Section:: *)
(*O(N) rank-2 anomalous dimensions*)


(* ::Subsection:: *)
(*Individual contribution from a single t-channel blocks to anomalous dimensions*)


(* ::Subsubsection:: *)
(*d=2 (AdS\:2083)*)


(* evaluate the general formula at dimensions of non-MFT (singlet) operators  *)
\[Gamma]tBlockList2d[{tnstart_,tnend_},\[CapitalDelta]ext_,\[Lambda]_,{n_,J_}] := \[Gamma]tBlockList2d[{tnstart,tnend},\[CapitalDelta]ext,\[Lambda],{n,J}] =
	With[
		{\[Gamma]tBlock = \[Gamma]tBlock2d[\[CapitalDelta]ext,{n,J},\[CapitalDelta]], poles = polesList2d[\[CapitalDelta]ext,\[Lambda]][[tnstart;;tnend]]},
		Parallelize[(\[Gamma]tBlock/.{\[CapitalDelta]->#}) & /@ poles]
	]
\[Gamma]tBlockList2d[nmax_,args___] := \[Gamma]tBlockList2d[{1,nmax},args]
\[Gamma]tBlockList2d[{tnstart_,tnend_},\[CapitalDelta]ext_,\[Lambda]_,J_] := \[Gamma]tBlockList2d[{tnstart,tnend},\[CapitalDelta]ext,\[Lambda],{0,J}]


(* ::Code:: *)
(*\[Gamma]tBlockList2d[200,7/4,10,0] // EchoTiming;*)
(*N[%,5]*)
(*(* the contributions grow, but we also need to weight them by OPE^2 coefficients *)*)
(*ListLogPlot[Abs@%]*)


(* ::Subsubsection:: *)
(*d=4 (AdS\:2085)*)


(* evaluate the general formula at dimensions of non-MFT (singlet) operators  *)
\[Gamma]tBlockList4d[{tnstart_,tnend_},\[CapitalDelta]ext_,\[Lambda]_,{n_,J_}] := \[Gamma]tBlockList4d[{tnstart,tnend},\[CapitalDelta]ext,\[Lambda],{n,J}] =
	With[
		{\[Gamma]tBlock = \[Gamma]tBlock4d[\[CapitalDelta]ext,{n,J},\[CapitalDelta]], poles = polesList4d[\[CapitalDelta]ext,\[Lambda]][[tnstart;;tnend]]},
		Parallelize[(\[Gamma]tBlock/.{\[CapitalDelta]->#}) & /@ poles]
	]
\[Gamma]tBlockList4d[nmax_,args___] := \[Gamma]tBlockList4d[{1,nmax},args]
\[Gamma]tBlockList4d[{tnstart_,tnend_},\[CapitalDelta]ext_,\[Lambda]_,J_] := \[Gamma]tBlockList4d[{tnstart,tnend},\[CapitalDelta]ext,\[Lambda],{0,J}]


(* ::Code:: *)
(*\[Gamma]tBlockList4d[50,7/3,10,0] // EchoTiming;*)
(*N[%,5]*)


(* ::Code:: *)
(*(* a little bit problematic sometimes, maybe differently simplified formula will be more robust *)*)
(*(* nonetheless, we can always choose more generic \[CapitalDelta]ext *)*)
(*\[Gamma]tBlock4d[\[CapitalDelta]ext,{0,J},\[CapitalDelta]p] // Simplify*)
(*\[Gamma]tBlock4d[\[CapitalDelta]ext,{n,J},\[CapitalDelta]p] /. {\[CapitalDelta]ext->9/2,n->0} // Simplify*)
(*\[Gamma]tBlock4d[9/2,{0,J},\[CapitalDelta]p] // Simplify*)
(*(\[Gamma]tBlock4d[\[CapitalDelta]ext,{0,J},\[CapitalDelta]p] // Simplify) /. \[CapitalDelta]ext->9/2 // Simplify*)


(* ::Subsection:: *)
(*Anomalous dimensions as weighted sum of t-channel contributions*)


(* ::Subsubsection:: *)
(*d=2  (AdS\:2083)*)


\[Gamma]List2d[{tnstart_,tnend_},\[CapitalDelta]ext_,\[Lambda]_,{n_,J_}] := \[Gamma]List2d[{tnstart,tnend},\[CapitalDelta]ext,\[Lambda],{n,J}] =
	With[
		{
			\[Gamma]tBlockList = \[Gamma]tBlockList2d[{tnstart,tnend},\[CapitalDelta]ext,\[Lambda],{n,J}],
			opeSqList = opeSqList2d[{tnstart,tnend},\[CapitalDelta]ext,\[Lambda]]
		},
		Parallelize[Inner[Times, \[Gamma]tBlockList, opeSqList, List]]
	]
\[Gamma]List2d[tnmax_,args___] := \[Gamma]List2d[{1,tnmax},args]
\[Gamma]List2d[{tnstart_,tnend_},\[CapitalDelta]ext_,\[Lambda]_,J_] := \[Gamma]List2d[{tnstart,tnend},\[CapitalDelta]ext,\[Lambda],{0,J}]


(* ::Code:: *)
(*AnomDimTable = Table[\[Gamma]List2d[20,7/4,10,{n,J}],{n,0,3},{J,0,3}] // EchoTiming;*)
(*Table[ListPlot[AnomDimTable[[n+1,J+1]], Filling->Axis,PlotRange->All,PlotLabel->StringTemplate["n=``,J=``"][n,J]],{n,0,3},{J,0,3}]*)
(*(* merge plots with the same n *)*)
(*Table[ListPlot[AnomDimTable[[n+1]], Filling->Axis,PlotLabel->StringTemplate["n=``"][n]],{n,0,3}]*)
(*(* merge plots with the same J *)*)
(*Table[ListLinePlot[AnomDimTable[[All,J+1]], Filling->Axis,PlotRange->{{0,10},All},PlotLabel->StringTemplate["J=``"][J]],{J,0,3}]*)


(* total anomalous dimensions of O(N) rank-2 {n,J} double-twist operators *)
\[Gamma]2d[{tnstart_,tnend_},\[CapitalDelta]ext_,\[Lambda]_,{n_,J_}] := \[Gamma]2d[{tnstart,tnend},\[CapitalDelta]ext,\[Lambda],{n,J}] =
	Plus@@\[Gamma]List2d[{tnstart,tnend},\[CapitalDelta]ext,\[Lambda],{n,J}]
\[Gamma]2d[tnmax_,args___] := \[Gamma]2d[{1,tnmax},args]
\[Gamma]2d[{tnstart_,tnend_},\[CapitalDelta]ext_,\[Lambda]_,J_] := \[Gamma]2d[{tnstart,tnend},\[CapitalDelta]ext,\[Lambda],{0,J}]


(* ::Code:: *)
(*(* dependence on tnend is not that relevant for J>0 *)*)
(*\[CapitalDelta]ext=E; \[Lambda]=E^7; nmax=3;*)
(*\[Gamma]Table         = Table[{n,J,\[Gamma]2d[50,\[CapitalDelta]ext,\[Lambda],{n,J}]},{n,0,nmax},{J,0,7}] // N // EchoTiming;*)
(*\[Gamma]TableSmaller  = Table[{n,J,\[Gamma]2d[30,\[CapitalDelta]ext,\[Lambda],{n,J}]},{n,0,nmax},{J,0,7}] // N // EchoTiming;*)
(*\[Gamma]TableSmallest = Table[{n,J,\[Gamma]2d[ 5,\[CapitalDelta]ext,\[Lambda],{n,J}]},{n,0,nmax},{J,0,7}] // N // EchoTiming;*)
(*(* differences most noticable for J=0 *)*)
(*(\[Gamma]Table-\[Gamma]TableSmaller)/Part[\[Gamma]Table,All,All,3]*)
(*(\[Gamma]Table-\[Gamma]TableSmallest)/Part[\[Gamma]Table,All,All,3]*)


(* ::Code:: *)
(*(* simple twist-spin plot *)*)
(*largeN=1/20;*)
(*Part[\[Gamma]Table,All,3]*)
(*twistTable = Apply[{#2,2\[CapitalDelta]ext+2#1+(1/largeN)#3}&,\[Gamma]Table,{2}]*)
(*MFTtwistTable = Apply[{#2,2\[CapitalDelta]ext+2#1}&,\[Gamma]Table,{2}]*)
(*twistPlot=ListPlot[twistTable];*)
(*twistLinePlot=ListLinePlot[twistTable,PlotStyle->Directive[Opacity[0.3],Thickness[0.001]],Filling->Table[n->2\[CapitalDelta]ext+2n-2,{n,1,20}],FillingStyle->Opacity[0.05]];*)
(*MFTtwistPlot=ListLinePlot[MFTtwistTable,PlotStyle->Directive[Dashing[0.01],GrayLevel[0.7],Thickness[0.001]]];*)
(*Show[MFTtwistPlot,twistLinePlot,twistPlot,ImageSize->Full]*)


(* ::Subsubsection::Closed:: *)
(*Simple testing of convergence properties*)


(* ::Code:: *)
(*n0J0List=\[Gamma]List2d[400,E,E^4,{0,0}];*)
(*ListLogLogPlot[Abs@%,AspectRatio->1]*)
(*a=399; b=390;*)
(*(Log@Abs@n0J0List[[a]]-Log@Abs@n0J0List[[b]])/(Log[a]-Log[b])*)


(* ::Code:: *)
(*n1J2List4d=\[Gamma]List4d[100,9/4,100,{1,2}];*)
(*ListLogLogPlot[Abs@%,AspectRatio->1]*)
(*a=99; b=90;*)
(*(Log@Abs@n1J2List4d[[a]]-Log@Abs@n1J2List4d[[b]])/(Log[a]-Log[b])*)


(* ::Code:: *)
(*n1J0List=\[Gamma]List2d[200,7/4,100,{1,0}];*)
(*ListLogLogPlot[Abs@%,AspectRatio->1]*)
(*a=199; b=150;*)
(*(Log@Abs@n1J0List[[a]]-Log@Abs@n1J0List[[b]])/(Log[a]-Log[b])*)


(* ::Code:: *)
(*n0J1List=\[Gamma]List2d[100,7/4,100,{0,1}];*)
(*ListLogLogPlot[Abs@%,AspectRatio->1]*)
(*a=99; b=80;*)
(*(Log@Abs@n0J1List[[a]]-Log@Abs@n0J1List[[b]])/(Log[a]-Log[b])*)


(* ::Code:: *)
(*n0J2List=\[Gamma]List2d[100,7/4,100,{0,2}];*)
(*ListLogLogPlot[Abs@%,AspectRatio->1]*)
(*a=99; b=80;*)
(*(Log@Abs@n0J2List[[a]]-Log@Abs@n0J2List[[b]])/(Log[a]-Log[b])*)


(* ::Code:: *)
(*n1J2List=\[Gamma]List2d[100,7/4,100,{1,2}];*)
(*ListLogLogPlot[Abs@%,AspectRatio->1]*)
(*a=99; b=80;*)
(*(Log@Abs@n1J2List[[a]]-Log@Abs@n1J2List[[b]])/(Log[a]-Log[b])*)


(* ::Code:: *)
(*n0J3List=\[Gamma]List2d[100,7/4,100,{0,3}];*)
(*ListLogLogPlot[Abs@%,AspectRatio->1]*)
(*a=99; b=80;*)
(*(Log@Abs@n0J3List[[a]]-Log@Abs@n0J3List[[b]])/(Log[a]-Log[b])*)


(* ::Subsubsection:: *)
(*d=4  (AdS\:2085)*)


\[Gamma]List4d[{tnstart_,tnend_},\[CapitalDelta]ext_,\[Lambda]_,{n_,J_}] := \[Gamma]List4d[{tnstart,tnend},\[CapitalDelta]ext,\[Lambda],{n,J}] =
	With[
		{
			\[Gamma]tBlockList = \[Gamma]tBlockList4d[{tnstart,tnend},\[CapitalDelta]ext,\[Lambda],{n,J}],
			opeSqList = opeSqList4d[{tnstart,tnend},\[CapitalDelta]ext,\[Lambda]]
		},
		Parallelize[Inner[Times, \[Gamma]tBlockList, opeSqList, List]]
	]
\[Gamma]List4d[tnmax_,args___] := \[Gamma]List4d[{1,tnmax},args]
\[Gamma]List4d[{tnstart_,tnend_},\[CapitalDelta]ext_,\[Lambda]_,J_] := \[Gamma]List4d[{tnstart,tnend},\[CapitalDelta]ext,\[Lambda],{0,J}]


(* total anomalous dimensions of O(N) rank-2 {n,J} double-twist operators *)
\[Gamma]4d[{tnstart_,tnend_},\[CapitalDelta]ext_,\[Lambda]_,{n_,J_}] := \[Gamma]4d[{tnstart,tnend},\[CapitalDelta]ext,\[Lambda],{n,J}] =
	Plus@@\[Gamma]List4d[{tnstart,tnend},\[CapitalDelta]ext,\[Lambda],{n,J}]
\[Gamma]4d[tnmax_,args___] := \[Gamma]4d[{1,tnmax},args]
\[Gamma]4d[{tnstart_,tnend_},\[CapitalDelta]ext_,\[Lambda]_,J_] := \[Gamma]4d[{tnstart,tnend},\[CapitalDelta]ext,\[Lambda],{0,J}]


(* ::Section::Closed:: *)
(*Epilog*)


End[]; (* End `Private` Context. *)
(* Protect[Evaluate[Context[] <> "*"]]; (* Protect all public symbols in the package. *) *)
EndPackage[];  (* End package Context. *)
