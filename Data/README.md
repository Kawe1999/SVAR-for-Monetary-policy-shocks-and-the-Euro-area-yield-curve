# Data README

All data is from ECB and Eurostat for the Euro Area and is already processed:

- OPR is the overnight policy rate and is spliced from ECB: EON.D.EONIA_TO.RATE and EST.B.EU000A2X2A25.WT
- y is hp-filtered industrial production from ECB: STS.M.I9.W.PROD.NS0020.4.000
- pi is HICP monthly data (annual rate of change) from Eurostat: prc_hicp_manr
- mo3-y30 is ECB zero-coupon yield curve data: 3month: YC.B.U2.EUR.4F.G_N_A.SV_C_YM.SR_3M , 1year: YC.B.U2.EUR.4F.G_N_A.SV_C_YM.IF_1Y , 2year: YC.B.U2.EUR.4F.G_N_A.SV_C_YM.SR_2Y , 5year: YC.B.U2.EUR.4F.G_N_A.SV_C_YM.SR_5Y , 7year: YC.B.U2.EUR.4F.G_N_A.SV_C_YM.SR_7Y , 10year: YC.B.U2.EUR.4F.G_N_A.SV_C_YM.IF_10Y , 20year: YC.B.U2.EUR.4F.G_N_A.SV_C_YM.SR_20Y , 30year YC.B.U2.EUR.4F.G_N_A.SV_C_YM.SR_30Y
- m3 is broad money supply from: BSI.M.U2.Y.V.M30.X.I.U2.2300.Z01.A
- dm3 is m3 first-difference.


