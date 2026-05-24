SetDirectory[FileNameJoin[{DirectoryName[$InputFileName], "..", ".."}]];

Get["MalmstenIntegralSequences/Kernel/MalmstenIntegralSequences.wl"];

tests = {
  VerificationTest[SechPowerIntegral[1], Pi/2, TestID -> "sech power n=1"],
  VerificationTest[SechPowerIntegral[2], 1, TestID -> "sech power n=2"],
  VerificationTest[SechPowerIntegral[3], Pi/4, TestID -> "sech power n=3"],

  VerificationTest[
    FullSimplify[
      MalmstenSequence[1] ==
        Pi Log[Gamma[3/4] Sqrt[2 Pi]/Gamma[1/4]]
    ],
    True,
    TestID -> "Malmsten n=1"
  ],

  VerificationTest[
    FullSimplify[
      MalmstenSequence[4] ==
        -2 EulerGamma/3 - 7 Zeta[3]/(3 Pi^2) + 2 Log[Pi/4]/3
    ],
    True,
    TestID -> "Malmsten n=4"
  ],

  VerificationTest[
    FixedTwoCoefficients[6],
    {0, 126, 56, 21, 6, 1},
    TestID -> "fixed-two coefficients"
  ],

  VerificationTest[
    StaircaseCoefficients[6],
    {42, 42, 28, 14, 5, 1},
    TestID -> "staircase coefficients"
  ],

  VerificationTest[
    NestedChiCoefficients[5, {2, 2, 2, 2, 2}],
    {0, 35, 15, 5, 1},
    TestID -> "nested chi coefficients"
  ],

  VerificationTest[
    NestedChiIntegrandPolynomial[4, {2, 2, 2, 2}, x],
    15 Sech[x] - 10 Sech[x]^2 - 4 Sech[x]^3 - Sech[x]^4,
    TestID -> "nested integrand polynomial"
  ],

  VerificationTest[
    NestedFSymbolicRHS[4, {2, 2, 2, 2}],
    10 \[ScriptCapitalF][2] + 4 \[ScriptCapitalF][3] + \[ScriptCapitalF][4],
    TestID -> "symbolic F RHS"
  ],

  VerificationTest[
    EvaluateFRHS[NestedFSymbolicRHS[4, {2, 2, 2, 2}]],
    NestedFRHS[4, {2, 2, 2, 2}],
    TestID -> "EvaluateFRHS matches NestedFRHS"
  ],

  VerificationTest[
    FSequence[3],
    ChiSequence[3],
    TestID -> "FSequence equals ChiSequence"
  ],

  VerificationTest[
    FreeQ[
      MalmstenSequence[4],
      StieltjesGamma[__] | Derivative[1, 0][Zeta][__] | Derivative[1][Zeta][__]
    ],
    True,
    TestID -> "Malmsten n=4 reduced"
  ],

  VerificationTest[
    HurwitzZetaDifference[1, a, b],
    PolyGamma[0, b] - PolyGamma[0, a],
    TestID -> "Hurwitz zeta pole difference"
  ],

  VerificationTest[
    HurwitzZetaPrimeDifference[0, a, b],
    LogGamma[a] - LogGamma[b],
    TestID -> "Hurwitz zeta prime zero difference"
  ]
};

report = TestReport[tests];
failed = Quiet[Check[report["TestsFailed"], Missing["Unavailable"]]];

Which[
  failed === {},
    Print["PASS: all ", Length[tests], " package checks passed."],
  ListQ[failed],
    Print["FAIL: ", Length[failed], " package check(s) failed."],
  True,
    Print["Test report generated. If the success rate is 100%, all package checks passed."]
];

report
