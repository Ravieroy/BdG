module BdgUtilities
    using LatticeUtilities
    using DelimitedFiles
    using PyCall
    using LinearAlgebra
    using CSV
    using Distributions
    using DataFrames
    using Random
    @pyimport numpy as np
    @pyimport random
    include("../src/generalutils.jl")
    using .GeneralUtils

    #Module exports
    export create_lattice_df
    export create_nnmat_df
    export create_t_matrix
    export create_bdg_ham
    export sort_evecs
    export check_norm
    export calc_delta
    export calc_avg_n
    export check_rel_tol
    export run_self_consistency_numpy
    export get_power_law_i_kagome
    export get_random_V
    export make_random_V_dict
    export get_scaled_V
    export get_corr_V_dict
    export get_effective_V
    export add_partial_disorder
    export create_partial_disorder_dict
    export add_sublattice_disorder
    export create_sublattice_disorder_dict
    export create_uncorrelated_disorder

    # alias to fermi function
    const f = fermi_fn

    """
        create_lattice_df(nnMap, unitCell, lattice, nSites; saveInFile=nothing, fileName=nothing)
    Returns a dataframe containing all the informtion relatd to lattice.
    """
    function create_lattice_df(nnMap, unitCell, lattice, nSites; saveInFile=nothing, fileName=nothing)
        df = DataFrame(siteIndex = Int[],
                    sitePosX = Float64[],
                    sitePosY = Float64[],
                    first = Int[],
                    second = Int[],
                    third = Int[],
                    fourth = Int[]
                    )
        for site in 1:nSites
            arr = zeros(ncol(df))
            loc = site_to_loc(site, unitCell, lattice)
            n1 = loc[1][1]
            n2 = loc[1][2]
            subLat = loc[2]
            pos = loc_to_pos([n1, n2], subLat, unitCell)
            arr[1] = site
            arr[2] = round(pos[1], digits=3)
            arr[3] = round(pos[2], digits=3)
            for i in 1:length(nnMap[1][2])
                arr[i+3] = nnMap[site][2][i]
            end
            push!(df, arr)
        end
        if saveInFile == true
            if fileName !== nothing
                CSV.write(fileName, df)
            else
                fileName = "../data/NN_MAP.csv"
                println("saving as $fileName")
                CSV.write(fileName, df)
            end
        end
        return df
    end

    """
        create_nnmat_df(N, lattice, bc; saveAsCsv=nothing, fileName=nothing)
    Returns the nearest neighbor map and dataframe of the given lattice.
    The dataframe can be saved as csv file.
    """
    function create_nnmat_df(N, lattice, bc; saveAsCsv=nothing, fileName=nothing)
        if lattice=="square"
            nSites = N ^ 2
            latticeVecs = [[1, 0.],[0., 1]]
            basisVecs = [[0.,0.]]
            square = UnitCell(latticeVecs, basisVecs)
            if bc == "pbc"
                periodic = [true, true]
            elseif bc == "obc"
                periodic = [false, false]
            else
                throw("Wrong boundary condition.(pbc || obc)")
            end
            L = [N, N]
            lat = Lattice(L, periodic)
            bondX = Bond(orbitals = (1,1), displacement = [1,0])
            bondY = Bond((1,1), [0,1])
            nbrTable = build_neighbor_table([bondX,bondY], square, lat)
            nbrTableMap = map_neighbor_table(nbrTable)
            df = create_lattice_df(nbrTableMap, square, lat, nSites; saveInFile=saveAsCsv, fileName=fileName)
            return square, lat, nbrTableMap, df

        elseif lattice=="kagome"
            nSites = 3 * N ^ 2
            latticeVecs = [[1.0,0.0], [1/2,√3/2]]
            basisVecs = [[0.0,0.0], [1/2,0.0], [1/4,√3/4]]
            kagome = UnitCell(latticeVecs, basisVecs)
            if bc == "pbc"
                periodic = [true, true]
            elseif bc == "obc"
                periodic = [false, false]
            else
                throw("Wrong boundary condition.(pbc || obc)")
            end
            L = [N, N]
            lat = Lattice(L, periodic)
            bond_1 = Bond(orbitals = (1,2), displacement = [0,0])
            bond_2 = Bond((1,3), [0,0])
            bond_3 = Bond((2,3), [0,0])
            bond_4 = Bond((2,1), [1,0])
            bond_5 = Bond((3,1), [0,1])
            bond_6 = Bond((3,2), [-1,1])
            nbrTable = build_neighbor_table([bond_1, bond_2, bond_3, bond_4, bond_5, bond_6], kagome, lat)
            nbrTableMap = map_neighbor_table(nbrTable)
            df = create_lattice_df(nbrTableMap, kagome, lat, nSites; saveInFile=saveAsCsv, fileName=fileName)
            return kagome, lat, nbrTableMap, df
        else
            throw("No implementation for $lattice(square || kagome)")
        end
    end


    """
        create_t_matrix(N, lattice, bc; saveAsText=tMatSave, fileName=tMatFileName)
    Return hopping matrix for the given `lattice` and its parameters.
    The t-matrix can be saved as text locally by setting `saveAsText=true`
    with `fileName=tMatFileName`.

    """
    function create_t_matrix(N, lattice, bc, nnMap; t=1, saveAsText=tMatSave, fileName=tMatFileName)
        nSites = length(nnMap)
        H = zeros(Float64, nSites, nSites)
        for site in 1: nSites
            for nbr in nnMap[site][2]
                H[site, nbr] = -t
            end
        end
        if saveAsText == true
            writedlm(fileName, H)
        end
        return H
    end

    """
        create_bdg_ham(deltaList, H, lattice, mu, N, nAvgOld, U, vList; nLayer=1)
    Returns the BdG Hamiltonian for given lattice.
    """
    function create_bdg_ham(deltaList, H, lattice, mu, N, nAvgOld, U, vList; nLayer=1)
        if lattice == "square"
            nSites = N ^ 2
        elseif lattice == "kagome"
            nSites = nLayer * 3 * N ^ 2
        else
            throw("No implementation for $lattice(square || kagome)")
        end

        for i in 1:nSites
            H[i, i] = -mu - (U/2) * nAvgOld[i] + vList[i]
        end

        # Making BdG Hamiltonian
        HBdG = zeros(Float64, 2 * nSites, 2 * nSites)
        for i in 1:nSites
            for j in 1:nSites
                HBdG[i, j] = H[i, j]
            end
        end

        for i in nSites+1: 2 * nSites
            for j in nSites+1: 2 * nSites
                HBdG[i, j] = -conj(H[i-nSites, j-nSites])
            end
        end

        # Add delta term
        for i in 1:nSites
            HBdG[i, i + nSites] = deltaList[i]
        end

        for i in nSites+1: 2 * nSites
            HBdG[i, i - nSites] = deltaList[i - nSites]
        end
        return HBdG
    end

    """
        sort_evecs(evector, nSites)
    Julia gives the eigenvectors in columns.
    i.e. First column is eigenvector for first eigenvalue and so on.
    evecs[:, n] gives the eigenvectors for nth eigenvalue.
    For BdG we want En>0 which starts from nSites+1 to 2*nSites.

    What we want is that all the rows contain eigenvectors.
    """
    function sort_evecs(evector, nSites)
        # evectors for En > 0 in the rows instead of colums as before.
        # So after the following line we will have un in 1:nSites column
        # and vn in nSites+1 : 2*nSites
        evecs = transpose(evector)[nSites + 1 : 2 * nSites, :]
        un = zeros(Float64, nSites, nSites)
        vn = zeros(Float64, nSites, nSites)

        #The following does a simple thing that it extracts the first
        # half of evecs as un and another half as vn.
        # There might be a way to directly extract the columns say,
        # from 1 to nSites.
        un = transpose(transpose(evecs)[1 : nSites, : ])
        vn = transpose(transpose(evecs)[nSites + 1 : 2 * nSites, : ])
        return un, vn
    end


    """
        check_norm(un, vn, nSites)
    Checks norm |un^2| + |vn^2| = 1 at every site for each eigenvalue, n
    """
    function check_norm(un, vn, nSites)
        normList = zeros(Float64, nSites)
        for n in 1:size(un)[1] # gives the number of rows in un, i.e. n
            val = dot(un[n, :], un[n, :]) + dot(vn[n, :], vn[n, :])
            normList[n] = val
        end
        return round.(normList, digits=5)
    end

    """
        calc_delta(U, un, vn, nSites)
    Calculates delta for the given un and vn
    """
    function calc_delta(U, un, vn, nSites, evals, T)
        deltaList = zeros(Float64, nSites)
        for i in 1:nSites
            delta = 0
            for n in 1:size(un)[1] # gives the number of rows in un, i.e. nth eigenvalue
                E = evals[n]
                # fn = f(E, T)
                fn = f(E, T=T)
                delta += U * (un[n, i] * conj(vn[n, i])) * (1 - 2 * fn)
                # delta += U * (un[n, i] * conj(vn[n, i]))
            end
            deltaList[i] = delta
        end
        return deltaList
    end

    """
        calc_delta(U, un, vn, nSites)
    Calculates N-average for the given vn
    """
    function calc_avg_n(un, vn, nSites, evals, T)
        nAvgList = zeros(Float64, nSites)
        for i in 1:nSites
            nAvg = 0
            for n in 1:size(vn)[1] # gives the number of rows in un, i.e. nth eigenvalue
                E = evals[n]
                # fn = f(E, T)
                fn = f(E, T=T)
                nAvg += ((un[n, i] * conj(un[n, i])) * fn + ((vn[n, i] * conj(vn[n, i])) * (1 - 2 * fn)))
                # nAvg += (vn[n, i] * conj(vn[n, i]))
            end
            nAvgList[i] = nAvg
        end
        return 2 * nAvgList
    end

    """
        run_self_consistency_numpy(deltaOld, mu, N, nAvgOld, nExp, lattice, seed, U, V; t = 1, tol=0.001, nLayer=1)
    Runs the self-consistent loop for BdG Hamiltonian.
    """
    function run_self_consistency_numpy(T, deltaOld, mu, N, nAvgOld, nExp, lattice, tMat, U, vList; tol=0.001, nLayer=1, maxCount=300)
        startTime = time()
        nSites = if lattice=="square" N^2 else nLayer * 3 * N ^ 2 end
        count = 0
        flag = false
        while flag == false
            count += 1
            H = create_bdg_ham(deltaOld, tMat, lattice, mu, N, nAvgOld, U, vList, nLayer=nLayer)
            (evals, evecs) = np.linalg.eigh(H)
            un, vn = sort_evecs(evecs, nSites)
            evalPostive = evals[nSites + 1 : end]
            deltaNew = calc_delta(U, un, vn, nSites, evalPostive, T)
            nAvgNew = calc_avg_n(un, vn, nSites, evalPostive, T)
            nAvg = sum(nAvgNew)/nSites
            mu += 0.3 * (nExp - nAvg)
            deltaFlag = check_rel_tol(deltaOld, deltaNew, tol=tol)
            nFlag = check_rel_tol(nAvgOld, nAvgNew, tol=tol)
            nAvgFlag = check_rel_tol([nExp], [nAvg], tol=1e-4)
            if deltaFlag == true && nFlag == true && nAvgFlag == true
                deltaFinal = deltaNew
                nAvgFinal = nAvgNew
                eGap = minimum(evals[nSites + 1: 2 * nSites])
                deltaGap = sum(deltaFinal)/nSites
                endTime = round(time() - startTime, digits = 2)
                isConverged = true
                return deltaFinal, nAvgFinal, eGap, deltaGap, evecs, evals, count, endTime, isConverged
                break
            else
                deltaOld = deltaNew
                nAvgOld = nAvgNew
                if count >= maxCount
                    deltaFinal = deltaNew
                    nAvgFinal = nAvgNew
                    eGap = minimum(evals[nSites + 1: 2 * nSites])
                    deltaGap = sum(deltaFinal)/nSites
                    endTime = round(time() - startTime, digits = 2)
                    isConverged = false
                    flag = true
                    return deltaFinal, nAvgFinal, eGap, deltaGap, evecs, evals, count, endTime, isConverged
                    break
                end
            end
        end #while loop end
    end # function end


    """
       get_power_law_i_kagome(alpha, df, N, seed, site)
    Returns random potential for a site from the prescription
    in Communications Physics, 5, 177 (2022)
    """
    function get_power_law_i_kagome(alpha, df, N, seed, site)
        ri = [df[site, :].sitePosX, df[site, :].sitePosY]
        res = 0
        np.random.seed(seed)
        for jx in 1:N/2
            for jy in 1:N/2
                phij = np.random.uniform(0, 2 * pi)
                qj = ((2 * pi * jx)/N, (2 * pi * jy)/N)
                factor = norm(qj, 2) ^ (-alpha/2)
                res += factor * cos(dot(qj, ri) + phij)
            end
        end
        return res/(3 * N ^ 2)
    end

    """
       get_random_V(alpha, df, N, seed)
    Returns an array with correlated random numbers for a given
    lattice(df),alpha, N and seed
    """
    function get_random_V(alpha, df, N, seed)
        nSites = nrow(df)
        vRandomList = zeros(Float64, nSites)
        for site in 1:nSites
            res = get_power_law_i_kagome(alpha, df, N, seed, site)
            vRandomList[site] = res
        end
        return vRandomList
    end

    """
        make_random_V_dict(seedList, alphaList, N, df)
    Returns a dictionary containing correlated random disorders
    (not normalized) for an array of seed and alpha.
    dict[seed][alpha]
    """
    function make_random_V_dict(seedList, alphaList, N, df)
        vDict = Dict()
        for seed in seedList
            vDictAlpha = Dict()
            for alpha in alphaList
                vList = get_random_V(alpha, df, N, seed)
                vDictAlpha[alpha] = vList
            end # alpha loop end
            vDict[seed] = vDictAlpha
        end # seed loop end
        return vDict
    end #function end

    """
        get_scaled_V(vList, V)
    Returns the scaled list of correlated random disorders for the given V
    """
    function get_scaled_V(vList, V)
        rMin = minimum(vList)
        rMax = maximum(vList)

        tMin = -V
        tMax = V

        vNorm = (vList .- rMin)./(rMax - rMin) .* (tMax - tMin) .+ tMin
        return vNorm
    end

    """
       get_corr_V_dict(seedList, alphaList, vList, N, df)
    Returns a dictionary containing correlated random disorders
    (normalized) for an array of seed, alpha and V.
    dict[seed][alpha][V]
    """
    function get_corr_V_dict(seedList, alphaList, vList, N, df)
        vDict = make_random_V_dict(seedList, alphaList, N, df)
        vCorrSeed = Dict()
        count = 0
        effVDict = Dict()
        for seed in seedList
            vCorrAlpha = Dict()
            for alpha in alphaList
                Vi = vDict[seed][alpha]
                vCorrV = Dict()
                for V in vList
                    count += 1
                    effV = get_effective_V(Vi, V)
                    #---- Storing ----------
                    effArr =  [seed, V, alpha, effV]
                    effVDict[count] = effArr
                    #------end storing------
                    vNorm = get_scaled_V(Vi, effV)
                    vCorr = vNorm .- sum(vNorm)/length(vNorm);
                    vCorrV[V] = vCorr
                end # V loop end
                vCorrAlpha[alpha] = vCorrV
            end # alpha loop end
            vCorrSeed[seed] = vCorrAlpha
        end # seed loop end
        effVname = "../data/effV"
        writedlm(effVname, sort(effVDict))
        open(effVname, "a") do file
            write(file, "count [seed, V, alpha, effV]")
        end
        return vCorrSeed
    end # function end

    """
        get_effective_V(vList, V)
    Returns the effective value for disorder for the given V.
    This function normalises and then fits the data with Normal distribution
    and returns the standard deviation.
    """
    function get_effective_V(vList, V)
        for val in 0:0.0001:15
            vNorm = get_scaled_V(vList, val)
            vCorr = vNorm .- sum(vNorm)/length(vNorm)
            fitParams = fit(Normal, vCorr)
            sigma = fitParams.σ
            if isapprox(sigma, V; rtol = 0.001)
                return val
                break
            end
        end
    end

    """
        add_partial_disorder(partialDisorderOn, V, seed, unitCell, lattice, N)
    Returns an array with random disorder on sites specified for the kagome lattice.
    """
    function add_partial_disorder(partialDisorderOn, V, seed, unitCell, lattice, N)
        nTrue = count(partialDisorderOn)
        nSites = 3 * N ^ 2
        np.random.seed(seed)
        val = np.random.uniform(-V, V, size=nTrue*N^2)
        vPartial = zeros(nSites)
        counter = 1
        for site in 1:nSites
            loc = site_to_loc(site, unitCell, lattice)
            subLat = loc[2]
            if partialDisorderOn[subLat] == true
                vPartial[site] = val[counter]
                counter += 1
            end # if block end
        end # site block end
        return vPartial
    end # function end


    """
        create_partial_disorder_dict(partialDisorderOn, vVals, seedList, unitCell, lattice, N)
    Returns a dictionary containing partial disorder for different V and seed.
    """
    function create_partial_disorder_dict(partialDisorderOn, vVals, seedList, unitCell, lattice, N)
        vPartialDict = Dict()
        for seed in seedList
            vPartialDictV = Dict()
            for V in vVals
                vPartialList = add_partial_disorder(partialDisorderOn, V, seed, unitCell, lattice, N)
                vPartialDictV[V] = vPartialList
            end # V block end
            vPartialDict[seed] = vPartialDictV
        end # seed block end
        return vPartialDict
    end

    """
        add_sublattice_disorder(subLatticeDisorder, V, seed, unitCell, lattice, N)
    Returns an array with disorder on random sites of sublattices with ratio defined
    in subLatticeDisorder.
    """
    function add_sublattice_disorder(subLatticeDisorder, V, seed, unitCell, lattice, N)
        nSites = 3 * N ^ 2
        TotalnumDisorderedSites = [0, 0, 0, Dict()]
        vPartial = zeros(nSites) # fill all the sites with V=0 initially

        if subLatticeDisorder.sites == "all"
            if subLatticeDisorder.input[1] == "ratio"
            numDisorderedSites = convert(Int, (subLatticeDisorder.input[2][2] * nSites))
            elseif subLatticeDisorder.input[1] == "number"
                numDisorderedSites = subLatticeDisorder.input[2][2]
            else
                throw("""
                        ERROR: valid options
                        input = ["ratio", [val1, val2]]
                        input = ["number", [n1, n2]]
                    """)
            end
            siteCounter = 1
            Random.seed!(seed)
            subLatSites = range(1, nSites, step=1) # creates a list of indices of sublattice
            siteIndices = sample(subLatSites, numDisorderedSites, replace = false) # chooses sites(random) to put disorder on
            # initialise/generate random potentials
            if subLatticeDisorder.disorderType == "random"
                Random.seed!(seed)
                vRandomList = rand(Uniform(-V, V), numDisorderedSites) # uniform random numbers between -V and V
            elseif subLatticeDisorder.disorderType == "constant"
                vRandomList = V * ones(numDisorderedSites)
            else
                throw("ERROR: choose either random or constant")
            end

            for site in siteIndices
                vPartial[site] = vRandomList[siteCounter]
                siteCounter += 1
                # get the type of lattice(A, B or C)
                loc = site_to_loc(site, unitCell, lattice) # get unit cell loc with site info
                if loc[2] == 1
                    TotalnumDisorderedSites[1] += 1
                elseif loc[2] == 2
                    TotalnumDisorderedSites[2] += 1
                else
                    TotalnumDisorderedSites[3] += 1
                end # if loc block end
            end #  for site block end
            TotalnumDisorderedSites[4] = sort(siteIndices) # sort and store the disordered site indices
        else # for disorder on specific lattice
            for subLat in subLatticeDisorder.sites
                # count the number of sites to put zero disorder on based on given ratio
                numDisorderedSites = convert(Int, (subLatticeDisorder.ratio[2] * (nSites/3)))
                siteCounter = 1
                Random.seed!(seed)
                subLatSites =range(subLat, nSites, step=3) # creates a list of indices of particular sublattice
                siteIndices = sample(subLatSites, numDisorderedSites, replace = false) # chooses sites(random) to put disorder on
                # initialise/generate random potentials
                if subLatticeDisorder.disorderType == "random"
                    Random.seed!(seed)
                    vRandomList = rand(Uniform(-V, V), numDisorderedSites) # uniform random numbers between -V and V
                elseif subLatticeDisorder.disorderType == "constant"
                    vRandomList = V * ones(numDisorderedSites)
                else
                    throw("ERROR: choose either random or constant")
                end
                for site in siteIndices
                    vPartial[site] = vRandomList[siteCounter]
                    TotalnumDisorderedSites[subLat] = siteCounter
                    siteCounter += 1
                end
                TotalnumDisorderedSites[4][subLat] = sort(siteIndices) # sort and store the disordered site indices
            end # subLat for loop block
        end # if block end
        return vPartial, TotalnumDisorderedSites
    end # function end
   """
       create_sublattice_disorder_dict(subLatticeDisorder, vVals, seed, unitCell, lattice, N)
   Returns a dict with an array of sublattice disorder for given seed list and V
   """
   function create_sublattice_disorder_dict(subLatticeDisorder, vVals, seedList, unitCell, lattice, N)
       vPartialDict = Dict()
       numDisorderedSitesDict = Dict()
       for seed in seedList
           vPartialDictV = Dict()
           numDisorderedSitesDictV = Dict()
           for V in vVals
               vPartialList, counter = add_sublattice_disorder(subLatticeDisorder, V, seed, unitCell, lattice, N)
               vPartialDictV[V] = vPartialList
               numDisorderedSitesDictV[V] = counter
           end # V block end
           vPartialDict[seed] = vPartialDictV
           numDisorderedSitesDict[seed] = numDisorderedSitesDictV
       end # seed block end
       return vPartialDict, numDisorderedSitesDict
   end


   """
        create_uncorrelated_disorder(N, V, seedList, vVals; pyRandom=true, nLayer=1)
   Returns a dictionary with uncorrelated disorders for the given seed and V values
   """
   function create_uncorrelated_disorder(N, lattice, seedList, vVals; pyRandom=true, nLayer=1)
       nSites = if lattice=="square" N ^ 2 else nLayer * 3 * N ^ 2 end
       vUncorrelatedDict = Dict()
       vList = zeros(nSites)
       for seed in seedList
           vUncorrelatedDictV = Dict()
           for V in vVals
               if pyRandom == true
                   np.random.seed(seed)
                   vList = np.random.uniform(-V, V, nSites)
                   vUncorrelatedDictV[V] = vList
               else
                   Random.seed!(seed)
                   vList = rand(Uniform(-V, V), nSites)
                   vUncorrelatedDictV[V] = vList
               end # pyRandom if block end
           end # V for loop end
           vUncorrelatedDict[seed] = vUncorrelatedDictV
       end # seed for loop end
       return vUncorrelatedDict
    end

end #module end
