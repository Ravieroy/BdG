module BdgUtilities
    using LatticeUtilities
    using DelimitedFiles
    using PyCall
    using LinearAlgebra
    @pyimport numpy as np

    include("../src/generalutils.jl")
    using .GeneralUtils

    #Module exports
    export create_t_matrix
    export create_nn_map 
    export create_bdg_ham
    export sort_evecs
    export check_norm
    export calc_delta
    export calc_avg_n
    export check_rel_tol
    export run_self_consistency_numpy
    """
        create_nn_map(N, lattice, bc; a=1, fileName=nnMapFileName, saveAsText=nnMapSave)
    Returns the nearest neighbors of given lattice as dict with keys as lattice
    and value as list of neighbors. Saves the neighbor table map as text 
    nbrTableMap=nbrTableMap, fileName=nnMapFileName)

    """
    function create_nn_map(N, lattice, bc; fileName=nnMapFileName, saveAsText=nnMapSave)
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
            if saveAsText == true
                save_nnmap_txt(nSites; nbrTableMap=nbrTableMap, fileName=fileName)
            end
            return nbrTableMap
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
            if saveAsText == true
                save_nnmap_txt(nSites; nbrTableMap=nbrTableMap, fileName=fileName)
            end
            return nbrTableMap
        else
            throw("No implementation for $lattice(square || kagome)")
        end 
    end

    """
        save_nnmap_txt(nSites; nbrTableMap=nbrTableMap, fileName=nnMapFileName)
    Saves the neighbor table map as text 
    nbrTableMap=nbrTableMap, fileName=nnMapFileName)

    """
    function save_nnmap_txt(nSites; nbrTableMap=nbrTableMap, fileName=nnMapFileName)
        open(fileName, "a") do io
            for i in 1:nSites
               writedlm(io, transpose(nbrTableMap[i][2]))
            end
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
        create_bdg_ham(deltaList, H, lattice, mu, N, nAvgOld, U, vList; t = 1)
    Returns the BdG Hamiltonian for given lattice.
    """
    function create_bdg_ham(deltaList, H, lattice, mu, N, nAvgOld, U, vList; t = 1)
    if lattice == "square"
        nSites = N ^ 2
    elseif lattice == "kagome"
        nSites = 3 * N ^ 2
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
    function calc_delta(U, un, vn, nSites)
        deltaList = zeros(Float64, nSites)
        for i in 1:nSites
            delta = 0
            for n in 1:size(un)[1] # gives the number of rows in un, i.e. nth eigenvalue
                delta += U * (un[n, i] * conj(vn[n, i]))
            end 
            deltaList[i] = delta
        end 
        return deltaList
    end 


    """
        calc_delta(U, un, vn, nSites)
    Calculates N-average for the given vn
    """
    function calc_avg_n(vn, nSites)
        nAvgList = zeros(Float64, nSites)
        for i in 1:nSites
            nAvg = 0
            for n in 1:size(vn)[1] # gives the number of rows in un, i.e. nth eigenvalue
                nAvg += (vn[n, i] * conj(vn[n, i]))
            end 
            nAvgList[i] = nAvg
        end 
        return 2 * nAvgList
    end

    """
        run_self_consistency_numpy(deltaOld, mu, N, nAvgOld, nExp, lattice, seed, U, V; t = 1, tol=0.001)
    Runs the self-consistent loop for BdG Hamiltonian.
    """
    function run_self_consistency_numpy(bc, deltaOld, mu, N, nAvgOld, nExp, lattice, seed, tMat, U, V; t = 1, tol=0.001)
        startTime = time()
        nSites = if lattice=="square" N^2 else 3*N^2 end
        count = 0
        np.random.seed(seed)
        vList = np.random.uniform(low=-V, high=V, size=(nSites,))
        flag = false
        while flag == false
            count += 1
            H = create_bdg_ham(deltaOld, tMat, lattice, mu, N, nAvgOld, U, vList; t = t)        
            (evals, evecs) = np.linalg.eigh(H)
            un, vn = sort_evecs(evecs, nSites)        
            deltaNew = calc_delta(U, un, vn, nSites)
            nAvgNew = calc_avg_n(vn, nSites)        
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
                println("Converged for V=$V, seed= $seed in $count iterations in $endTime s = $(round(endTime/60, digits=2)) mins")
                return deltaFinal, nAvgFinal, eGap, deltaGap, evecs, evals, count, endTime
                break           
            else
                deltaOld = deltaNew
                nAvgOld = nAvgNew
                flag = false
            end
        end #while loop end 
    end # function end 
end #module end  
