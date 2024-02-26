N = 12
lattice = "kagome"
nSites = 3 * N ^ 2
a = 1
t = 1
bc = "pbc"
tMatSave = false
tMatFileName = "tmat_$lattice$N"
nExp = 0.875
tol = 0.001
nAvgOld = ones(Float64, nSites)
deltaOld = ones(Float64, nSites)
mu = -0.92
U = 1.5 * t
seedList = [1, 3]
alphaList = [0, 1  ]
vVals = [0.1, 1]
flag = false
dataSetFolder = "../data/"
fileFormat = "pkl"
saveInFolder = "../results"
store = Dict([
    ("store_delta_op", [true, "deltaDict", 1]),
    ("store_avg_n", [true, "nAvgDict", 2]),
    ("store_egap", [true, "eGapDict", 3]),
    ("store_delta_gap", [true, "deltaGapDict", 4]),
    ("store_evectors", [true, "evectors", 5]),
    ("store_evalues", [true, "evalues", 6]),
])

