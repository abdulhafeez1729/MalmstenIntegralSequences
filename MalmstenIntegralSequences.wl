(* ::Package:: *)

BeginPackage["MalmstenIntegralSequences`"];

SignedStirlingP::usage =
  "SignedStirlingP[k, m, z] gives the signed generalized Stirling polynomial \
P_k(m,z) used in the integral-sequence formulae.";

SechPowerIntegral::usage =
  "SechPowerIntegral[n] gives Integrate[Sech[x]^n,{x,0,Infinity}]. \
SechPowerIntegral[n,b] gives the corresponding scaled value for b > 0.";

MalmstenSequence::usage =
  "MalmstenSequence[n,a,b] gives the reduced closed form for \
Integrate[Log[a x] Sech[b x]^n,{x,0,Infinity}] in the positive-parameter \
range a > 0, b > 0. MalmstenSequence[n] uses a=b=1.";

HurwitzZetaDifference::usage =
  "HurwitzZetaDifference[s,a,b] gives Zeta[s,a]-Zeta[s,b], with removable \
singularities simplified using Lemma 3-style identities.";

HurwitzZetaPrimeDifference::usage =
  "HurwitzZetaPrimeDifference[s,a,b] gives D[Zeta[s,a]-Zeta[s,b],s], with \
removable singularities simplified.  At s=1 it returns the finite part obtained \
from the Laurent expansion of the Hurwitz zeta function.";

HurwitzHalfDifference::usage =
  "HurwitzHalfDifference[s,x] gives Zeta[s,x]-Zeta[s,x+1/2].";

HurwitzHalfPrimeDifference::usage =
  "HurwitzHalfPrimeDifference[s,x] gives the s-derivative of \
Zeta[s,x]-Zeta[s,x+1/2].";

ChiSequence::usage =
  "ChiSequence[n] gives the reduced closed form for chi_n, where chi_n = \
Integrate[(Sech[x] - Sech[x]^n)/x^2, {x,0,Infinity}].";

FSequence::usage =
  "FSequence[j] gives the reduced closed-form term script F_j used in the \
nested-sum identities. It uses the same reduction as ChiSequence[j] and \
NestedFRHS[N,lows].";

LambdaSequence::usage =
  "LambdaSequence[n] gives lambda_n, using n lambda_n = chi_n + lambda_1.";

DeltaSequence::usage =
  "DeltaSequence[n] gives delta_n, using delta_n = chi_(n+1) - chi_n.";

NestedChiCoefficients::usage =
  "NestedChiCoefficients[N, {l1,...,lN}] returns the multiplicities \
{c1,...,cN} of chi_1,...,chi_N in the nested sum from equation (56).";

NestedFRHS::usage =
  "NestedFRHS[N, lows] returns the reduced closed-form right-hand side \
Sum[c_j FSequence[j], {j,1,N}] associated with the nested sum from equation \
(56).";

ReduceQuarterZetaTerms::usage =
  "ReduceQuarterZetaTerms[expr] reduces Hurwitz-zeta derivative terms at \
quarter arguments using shift recurrences and the complementary relation from \
Lemma 3, equation (22), of the paper.";

NestedFRHSReduced::usage =
  "NestedFRHSReduced[N, lows] is kept as an alias for NestedFRHS[N, lows].";

NestedFSymbolicRHS::usage =
  "NestedFSymbolicRHS[N, lows] returns the symbolic right-hand side \
Sum[c_j \[ScriptCapitalF][j], {j,1,N}], without expanding the closed forms. \
NestedFSymbolicRHS[N, lows, f] uses f as the symbolic head.";

EvaluateFRHS::usage =
  "EvaluateFRHS[expr] replaces symbolic script F terms in expr by the reduced \
values FSequence[j], then reduces the resulting expression. EvaluateFRHS[expr, \
f] uses f as the symbolic head.";

NestedChiIntegrandPolynomial::usage =
  "NestedChiIntegrandPolynomial[N, lows, x] returns the polynomial p(Sech[x]) \
such that the nested sum equals Integrate[p(Sech[x])/x^2, {x,0,Infinity}].";

FixedTwoCoefficients::usage =
  "FixedTwoCoefficients[N] gives the coefficient vector for l1=...=lN=2.";

StaircaseCoefficients::usage =
  "StaircaseCoefficients[N] gives the coefficient vector for \
{l1,...,lN}={1,2,...,N}.";

Begin["`Private`"];

ClearAll[zetaPrime];
zetaPrime[s_, a_] := Derivative[1, 0][Zeta][s, a];

(* These difference helpers prevent Mathematica from forming expressions such as
   ComplexInfinity - ComplexInfinity at removable singularities.  The s=1 rule
   for HurwitzZetaPrimeDifference is the finite part from the Laurent expansion
   of the Hurwitz zeta function, not a separately stated theorem in the paper. *)
ClearAll[HurwitzZetaDifference, HurwitzZetaPrimeDifference];
HurwitzZetaDifference[1, a_, b_] := PolyGamma[0, b] - PolyGamma[0, a];
HurwitzZetaDifference[s_, a_, b_] := Zeta[s, a] - Zeta[s, b];

HurwitzZetaPrimeDifference[1, 1/4, 3/4] := quarterPoleStieltjesDifference;
HurwitzZetaPrimeDifference[1, 3/4, 1/4] := -quarterPoleStieltjesDifference;
HurwitzZetaPrimeDifference[1, a_, b_] :=
  ReduceQuarterZetaTerms[-StieltjesGamma[1, a] + StieltjesGamma[1, b]];
HurwitzZetaPrimeDifference[0, a_, b_] := LogGamma[a] - LogGamma[b];
HurwitzZetaPrimeDifference[s_, a_, b_] := zetaPrime[s, a] - zetaPrime[s, b];

ClearAll[HurwitzHalfDifference, HurwitzHalfPrimeDifference];
HurwitzHalfDifference[s_, x_] := HurwitzZetaDifference[s, x, x + 1/2];
HurwitzHalfPrimeDifference[s_, x_] := HurwitzZetaPrimeDifference[s, x, x + 1/2];

ClearAll[cleanSimplifyArtifacts, shiftQuarterTerms, reduceHalfIntegerTerms,
  quarterPoleStieltjesDifference, reduceQuarterPoleTerms, reflectionTerm22,
  reduceReflectionPairs, ReduceQuarterZetaTerms];

cleanSimplifyArtifacts[expr_] :=
  expr //. {
    HoldPattern[Simplify[x_, Assumptions -> True]] :> x,
    HoldPattern[FullSimplify[x_, Assumptions -> True]] :> x
  };

shiftQuarterTerms[expr_] :=
  expr //. {
    StieltjesGamma[1, r_Rational] /; r > 1 :>
      Module[{n = Floor[r], a = r - Floor[r]},
        If[a == 0, a = 1; n--];
        StieltjesGamma[1, a] - Sum[Log[a + k]/(a + k), {k, 0, n - 1}]
      ],
    Derivative[1, 0][Zeta][s_Integer, r_Rational] /; r > 1 :>
      Module[{n = Floor[r], a = r - Floor[r]},
        If[a == 0, a = 1; n--];
        zetaPrime[s, a] + Sum[(a + k)^(-s) Log[a + k], {k, 0, n - 1}]
      ],
    Zeta[s_Integer, r_Rational] /; r > 1 :>
      Module[{n = Floor[r], a = r - Floor[r]},
        If[a == 0, a = 1; n--];
        Zeta[s, a] - Sum[(a + k)^(-s), {k, 0, n - 1}]
      ]
  };

reduceHalfIntegerTerms[expr_] :=
  expr //. {
    Zeta[s_Integer, 1] :> Zeta[s],
    Zeta[s_Integer, 1/2] :> (2^s - 1) Zeta[s],
    Derivative[1, 0][Zeta][s_Integer, 1] :> Derivative[1][Zeta][s],
    Derivative[1, 0][Zeta][s_Integer, 1/2] :>
      (2^s - 1) Derivative[1][Zeta][s] + 2^s Log[2] Zeta[s],
    Derivative[1][Zeta][0] :> -Log[2 Pi]/2,
    Derivative[1][Zeta][-2 r_Integer?Positive] :>
      (-1)^r Factorial[2 r] Zeta[2 r + 1]/(2 (2 Pi)^(2 r)),
    StieltjesGamma[1, 1] :> StieltjesGamma[1],
    StieltjesGamma[1, 1/2] :>
      StieltjesGamma[1] - 2 EulerGamma Log[2] - Log[2]^2
  };

quarterPoleStieltjesDifference :=
  Pi (EulerGamma + 4 Log[2] + 3 Log[Pi] - 4 Log[Gamma[1/4]]);

(* The reduction above handles zeta-prime pairs at nonpositive integers.
   This handles the corresponding pole finite part at s=1:
   -StieltjesGamma[1,1/4] + StieltjesGamma[1,3/4]. *)
reduceQuarterPoleTerms[expr_] :=
  Module[{e = Expand[expr], a, b, ca, cb},
    a = StieltjesGamma[1, 1/4];
    b = StieltjesGamma[1, 3/4];
    ca = FullSimplify[Coefficient[e, a]];
    cb = FullSimplify[Coefficient[e, b]];
    If[ca =!= 0 && FullSimplify[cb == -ca],
      e = FullSimplify[e - ca a - cb b + cb quarterPoleStieltjesDifference]
    ];
    e
  ];

reflectionTerm22[n_Integer?NonNegative, z_] :=
  Pi I BernoulliB[n + 1, z]/(n + 1)
  + Factorial[n] Exp[-Pi I n/2] PolyLog[n + 1, Exp[2 Pi I z]]/(2 Pi)^n;

reduceReflectionPairs[expr_] :=
  Module[{e = Expand[expr], ns, a, b, ca, cb, sign, rhs},
    ns = DeleteDuplicates @ Cases[
      e,
      Derivative[1, 0][Zeta][-n_Integer?NonNegative, 1/4 | 3/4] :> n,
      {0, Infinity}
    ];
    Do[
      a = zetaPrime[-n, 1/4];
      b = zetaPrime[-n, 3/4];
      ca = FullSimplify[Coefficient[e, a]];
      cb = FullSimplify[Coefficient[e, b]];
      sign = (-1)^n;
      If[ca =!= 0 && FullSimplify[cb == sign ca],
        rhs = reflectionTerm22[n, 1/4];
        e = FullSimplify[e - ca (a + sign b) + ca rhs]
      ],
      {n, ns}
    ];
    e
  ];

ReduceQuarterZetaTerms[expr_] :=
  Module[{e, atoms},
    e = cleanSimplifyArtifacts @ FullSimplify[
      cleanSimplifyArtifacts @ FunctionExpand[
        reduceHalfIntegerTerms @ reduceQuarterPoleTerms @ reduceReflectionPairs[
          reduceHalfIntegerTerms @ reduceQuarterPoleTerms @ shiftQuarterTerms[expr]
        ]
      ]
    ];
    atoms = DeleteDuplicates @ Cases[
      e,
      Catalan | Log[_] | Zeta[_] | PolyLog[_, _] | StieltjesGamma[__] |
        Derivative[1, 0][Zeta][_, _] | Derivative[1][Zeta][_],
      {0, Infinity}
    ];
    cleanSimplifyArtifacts @ FullSimplify[
      reduceHalfIntegerTerms @ Collect[e, atoms, FullSimplify]
    ]
  ];

ClearAll[SignedStirlingP];
SignedStirlingP[k_Integer?NonNegative, m_Integer?NonNegative, z_] /; k <= m :=
  Sum[
    (-1)^(i - k) z^i Binomial[i + m - k, m - k]
      Abs[StirlingS1[m + 1, i + m - k + 1]],
    {i, 0, k}
  ];

ClearAll[SechPowerIntegral];
SechPowerIntegral[n_Integer?Positive] :=
  SechPowerIntegral[n] =
    FunctionExpand[2^(n - 2) Gamma[n/2]^2/Factorial[n - 1]];

SechPowerIntegral[n_Integer?Positive, b_] :=
  SechPowerIntegral[n]/b;

ClearAll[malmstenSequenceRaw, MalmstenSequence];
MalmstenSequence[n_Integer?Positive] := MalmstenSequence[n, 1, 1];

malmstenSequenceRaw[1, a_, b_] := malmstenSequenceRaw[1, a, b] =
  FunctionExpand[
    Pi/b Log[Gamma[3/4] Sqrt[2 Pi] Sqrt[a/b]/Gamma[1/4]]
  ];

malmstenSequenceRaw[n_Integer?Positive, a_, b_] /; n >= 2 :=
  malmstenSequenceRaw[n, a, b] =
    Expand @ FunctionExpand[
      SechPowerIntegral[n, b] Log[a/b]
      + (2^(2 n - 1)/(b Factorial[n - 1]))
        Sum[
          (-1/2)^k SignedStirlingP[k - 2, n - 1, n/2]
            (
              HurwitzZetaPrimeDifference[k - n, n/4, (n + 2)/4]
              - (EulerGamma + Log[4])
                HurwitzZetaDifference[k - n, n/4, (n + 2)/4]
            ),
          {k, 2, n + 1}
        ]
    ];

MalmstenSequence[n_Integer?Positive, a_, b_] :=
  MalmstenSequence[n, a, b] =
    ReduceQuarterZetaTerms[malmstenSequenceRaw[n, a, b]];

ClearAll[chiSequenceRaw, ChiSequence];
chiSequenceRaw[1] = 0;
chiSequenceRaw[n_Integer?Positive] := chiSequenceRaw[n] =
  Expand @ FunctionExpand[
    -4 Catalan/Pi
    + (2^(2 n - 3) n^2/Factorial[n - 1])
      Sum[
        (-1/2)^k SignedStirlingP[k, n - 1, n/2]
          HurwitzZetaPrimeDifference[k - n + 2, n/4, (n + 2)/4],
        {k, 0, n - 1}
      ]
    - (2^(2 n + 1)/Factorial[n - 1])
      Sum[
        (-1/2)^k SignedStirlingP[k, n + 1, (n + 2)/2]
          HurwitzZetaPrimeDifference[k - n, (n + 2)/4, (n + 4)/4],
        {k, 0, n + 1}
      ]
  ];

ChiSequence[n_Integer?Positive] := ChiSequence[n] =
  ReduceQuarterZetaTerms[chiSequenceRaw[n]];

ClearAll[FSequence];
FSequence[n_Integer?Positive] := FSequence[n] =
  ReduceQuarterZetaTerms[ChiSequence[n]];

ClearAll[LambdaSequence];
LambdaSequence[1] = 4 Catalan/Pi;
LambdaSequence[n_Integer?Positive] := LambdaSequence[n] =
  ReduceQuarterZetaTerms @ Expand @ FunctionExpand[
    (ChiSequence[n] + LambdaSequence[1])/n
  ];

ClearAll[DeltaSequence];
DeltaSequence[n_Integer?Positive] := DeltaSequence[n] =
  ReduceQuarterZetaTerms @ Expand @ FunctionExpand[
    ChiSequence[n + 1] - ChiSequence[n]
  ];

ClearAll[NestedChiCoefficients];
NestedChiCoefficients::bounds =
  "Expected a list of exactly `1` positive integer lower bounds.";

NestedChiCoefficients[n_Integer?Positive, lows : {__Integer?Positive}] /;
    Length[lows] == n :=
  Module[{counts = ConstantArray[0, n], rec},
    rec[1, upper_Integer] := Do[
      If[1 <= k <= n, counts[[k]]++],
      {k, lows[[1]], upper}
    ];
    rec[level_Integer, upper_Integer] := Do[
      rec[level - 1, k],
      {k, lows[[level]], upper}
    ];

    If[n == 1,
      Do[If[1 <= k <= n, counts[[k]]++], {k, lows[[1]], n}],
      Do[rec[n - 1, k], {k, lows[[n]], n}]
    ];
    counts
  ];

NestedChiCoefficients[n_Integer?Positive, lows_] := (
  Message[NestedChiCoefficients::bounds, n];
  $Failed
);

ClearAll[nestedFRHSRaw];
nestedFRHSRaw[n_Integer?Positive, lows_List] :=
  Module[{c = NestedChiCoefficients[n, lows]},
    If[c === $Failed, $Failed, Total[MapThread[#1 FSequence[#2] &, {c, Range[n]}]]]
  ];

ClearAll[NestedFRHS];
NestedFRHS[n_Integer?Positive, lows_List] :=
  Module[{rhs = nestedFRHSRaw[n, lows]},
    If[rhs === $Failed, $Failed, ReduceQuarterZetaTerms[rhs]]
  ];

ClearAll[NestedFRHSReduced];
NestedFRHSReduced[n_Integer?Positive, lows_List] :=
  NestedFRHS[n, lows];

ClearAll[NestedFSymbolicRHS];
NestedFSymbolicRHS[n_Integer?Positive, lows_List, f_: Global`\[ScriptCapitalF]] :=
  Module[{c = NestedChiCoefficients[n, lows]},
    If[c === $Failed, $Failed, Total[MapThread[#1 f[#2] &, {c, Range[n]}]]]
  ];

ClearAll[replaceFSymbols, EvaluateFRHS];
(* EvaluateFRHS is deliberately two-step: substitute symbolic F_j terms, then
   apply the same reduction used by NestedFRHS. *)
replaceFSymbols[expr_, f_] :=
  Module[{h = f},
    expr /. {
      h[j_Integer?Positive] :> FSequence[j],
      Subscript[h, j_Integer?Positive] :> FSequence[j]
    }
  ];

EvaluateFRHS[expr_, f_: Global`\[ScriptCapitalF]] :=
  ReduceQuarterZetaTerms[replaceFSymbols[expr, f]];

ClearAll[NestedChiIntegrandPolynomial];
NestedChiIntegrandPolynomial[n_Integer?Positive, lows_List, x_: Global`x] :=
  Module[{c = NestedChiCoefficients[n, lows]},
    If[c === $Failed,
      $Failed,
      Total[c] Sech[x] - Total[MapThread[#1 Sech[x]^#2 &, {c, Range[n]}]]
    ]
  ];

ClearAll[FixedTwoCoefficients];
FixedTwoCoefficients[n_Integer?Positive] :=
  If[n == 1,
    {0},
    Prepend[Table[Binomial[2 n - j - 1, n - j], {j, 2, n}], 0]
  ];

ClearAll[StaircaseCoefficients];
StaircaseCoefficients[n_Integer?Positive] :=
  NestedChiCoefficients[n, Range[n]];

End[];

EndPackage[];
