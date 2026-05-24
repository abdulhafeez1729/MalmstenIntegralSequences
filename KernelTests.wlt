VerificationTest[
  Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "Kernel", "MalmstenIntegralSequences.wl"}]];
  True,
  True,
  TestID -> "Load package"
]

VerificationTest[
  FixedTwoCoefficients[6],
  {0, 126, 56, 21, 6, 1},
  TestID -> "Fixed lower bounds at two coefficients"
]

VerificationTest[
  StaircaseCoefficients[6],
  {42, 42, 28, 14, 5, 1},
  TestID -> "Staircase coefficients"
]

VerificationTest[
  NestedChiCoefficients[5, {2, 2, 2, 2, 2}],
  {0, 35, 15, 5, 1},
  TestID -> "Nested coefficients for all lower bounds two"
]

VerificationTest[
  NestedChiIntegrandPolynomial[4, {2, 2, 2, 2}, x],
  15 Sech[x] - 10 Sech[x]^2 - 4 Sech[x]^3 - Sech[x]^4,
  TestID -> "Integrand polynomial for N=4 all lower bounds two"
]

VerificationTest[
  SechPowerIntegral[1],
  Pi/2,
  TestID -> "Integral of sech x"
]

VerificationTest[
  SechPowerIntegral[2],
  1,
  TestID -> "Integral of sech^2 x"
]

VerificationTest[
  SechPowerIntegral[3],
  Pi/4,
  TestID -> "Integral of sech^3 x"
]

VerificationTest[
  FullSimplify[
    MalmstenSequence[1] ==
      Pi Log[Gamma[3/4] Sqrt[2 Pi]/Gamma[1/4]]
  ],
  True,
  TestID -> "Malmsten sequence n=1"
]

VerificationTest[
  HurwitzZetaDifference[1, a, b],
  PolyGamma[0, b] - PolyGamma[0, a],
  TestID -> "Hurwitz zeta difference at pole"
]

VerificationTest[
  HurwitzZetaPrimeDifference[1, a, b],
  -StieltjesGamma[1, a] + StieltjesGamma[1, b],
  TestID -> "Hurwitz zeta derivative difference at pole"
]

VerificationTest[
  FreeQ[HurwitzZetaPrimeDifference[1, 1/4, 3/4], StieltjesGamma[__]],
  True,
  TestID -> "Hurwitz zeta derivative quarter pole difference"
]

VerificationTest[
  ReduceQuarterZetaTerms[
    a (-StieltjesGamma[1, 1/4] + StieltjesGamma[1, 3/4])
  ],
  a Pi (EulerGamma + 4 Log[2] + 3 Log[Pi] - 4 Log[Gamma[1/4]]),
  TestID -> "Generic reduction of quarter pole finite part"
]

VerificationTest[
  HurwitzZetaPrimeDifference[0, a, b],
  LogGamma[a] - LogGamma[b],
  TestID -> "Hurwitz zeta derivative difference at zero"
]

VerificationTest[
  FullSimplify[
    MalmstenSequence[4] ==
      -2 EulerGamma/3 - 7 Zeta[3]/(3 Pi^2) + 2 Log[Pi/4]/3
  ],
  True,
  TestID -> "Malmsten sequence n=4"
]

VerificationTest[
  FreeQ[
    MalmstenSequence[4],
    StieltjesGamma[__] | Derivative[1, 0][Zeta][__] | Derivative[1][Zeta][__]
  ],
  True,
  TestID -> "MalmstenSequence returns reduced half-integer form"
]

VerificationTest[
  NestedChiIntegrandPolynomial[4, {1, 2, 3, 4}, x],
  9 Sech[x] - 5 Sech[x]^2 - 3 Sech[x]^3 - Sech[x]^4,
  TestID -> "Integrand polynomial for N=4 staircase"
]

VerificationTest[
  NestedFSymbolicRHS[4, {2, 2, 2, 2}, f],
  10 f[2] + 4 f[3] + f[4],
  TestID -> "Symbolic F RHS for N=4 all lower bounds two"
]

VerificationTest[
  NestedFSymbolicRHS[4, {1, 2, 3, 4}, f],
  5 f[1] + 5 f[2] + 3 f[3] + f[4],
  TestID -> "Symbolic F RHS for N=4 staircase"
]

VerificationTest[
  FreeQ[ChiSequence[3], StieltjesGamma[__] | Derivative[1, 0][Zeta][__]],
  True,
  TestID -> "ChiSequence returns reduced quarter form"
]

VerificationTest[
  FSequence[3],
  ChiSequence[3],
  TestID -> "FSequence uses reduced ChiSequence value"
]

VerificationTest[
  EvaluateFRHS[3 Global`\[ScriptCapitalF][2] + Global`\[ScriptCapitalF][3]],
  ReduceQuarterZetaTerms[3 FSequence[2] + FSequence[3]],
  TestID -> "EvaluateFRHS substitutes and reduces script F terms"
]
