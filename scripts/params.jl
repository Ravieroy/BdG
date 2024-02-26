N                       =  24
nnb                     =  4
lattice                 =  "square"
nLayer                  =  1
nSites                  =  nLayer * N ^ 2
a                       =  1
t                       =  1
bc                      =  "pbc"
tMatSave                =  false
tMatFileName            =  "tmat_$lattice$N"
nExp                    =  0.875
tol                     =  0.001
nAvgOld                 =  ones(Float64, nSites)
deltaOld                =  ones(Float64, nSites)
maxCount                =  200
mu                      =  -0.92
U                       =  1.5  * t
tempList                =  [0]
seedList                =  [1, 3, 5, 7, 9]
vVals                   =  [0.1, 0.5, 1.25, 2, 3]
alphaList               =  [0]
correlated              =  false
nnHam                   =  (type =  "!prebuilt", fileName =  "../data/ham_$N")
partial                 =  false
sublattice              =  false
uncorrelated            =  true
continuedCalc           = (state = false, variable = "V", run=2)
# partialDisorderOn       =  [true, true, true] # [siteA, siteB, siteC]
# subLatticeDisorder    =  (sites =  "all", input =  ["ratio", [9/10, 1/10]],  disorderType="random")
# subLatticeDisorder      =  (sites =  "all", input =  ["number", [270, 30]], disorderType="random")
dataSetFolder           =  "../data/"
logFolder               =  "../logs/"
fileFormat              =  "pkl"
saveInFolder            =  "../results"
store                   =  Dict([
    ("store_delta_op", [true, "deltaDict", 1]),
    ("store_avg_n", [true, "nAvgDict", 2]),
    ("store_egap", [true, "eGapDict", 3]),
    ("store_delta_gap", [true, "deltaGapDict", 4]),
    ("store_evectors", [true, "evectors", 5]),
    ("store_evalues", [true, "evalues", 6]),
])

