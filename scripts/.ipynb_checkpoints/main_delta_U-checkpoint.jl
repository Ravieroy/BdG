# Importing required libraries
using Distributed
@everywhere begin
    using LatticeUtilities
    using DelimitedFiles
    using LinearAlgebra
    using Dates
    using PyCall
    @everywhere(@pyimport numpy as np)
    @everywhere(@pyimport random)
    # Importing personal modules
    include("../src/bdg_utilities.jl")
    include("../src/generalutils.jl")
    include("params.jl")
end # module everywhere end

# Loading personal modules
@everywhere using .BdgUtilities
@everywhere using .GeneralUtils

@everywhere function main(U)
    #------Logging part-----------
    logFileName = string(logFolder, "info","_U_$U", "_id_$(myid())", ".log")
    dateFormat = "yyyy-mm-dd HH:MM:SS"
    startDate = Dates.format(now(), dateFormat)
    message1 = "\nBdG Calculation for N=$N for U=$U started : $startDate\n"
    message2 = string("number of processes : $(nprocs())\n", "number of workers : $(nworkers())\n" )
    message = string(message1, message2, "\n")
    write_to_file(message, logFileName)
    #-----------Logging part end ------

    deltaDictSeed = Dict()
    for seed in seedList
        deltaDictα = Dict()
        for alpha in alphaList
            deltaDictV = Dict()
            for V in vVals
                nAvgOld = ones(Float64, nSites)
                deltaOld = ones(Float64, nSites)
                mu = -0.9
                if V == 0
                    vList = zeros(nSites)
                else 
                    vList = vDict[seed][alpha][V]
                end 
                #  Logging part
                timeNow = Dates.format(now(), "HH:MM")
                message = "($timeNow) Started for seed=$seed, V=$V, U=$U\n"
                write_to_file(message, logFileName)
                deltaFinal, _, _, _, _, _, count, endTime, isConverged =
                                    run_self_consistency_numpy(deltaOld,
                                                                mu,
                                                                N,
                                                                nAvgOld,
                                                                nExp,
                                                                lattice,
                                                                seed,
                                                                tMat,
                                                                U,
                                                                vList,
                                                                tol=tol,
                                                                nLayer=nLayer,
                                                               maxCount=maxCount)
                if isConverged == true 
                    timeNow = Dates.format(now(), "HH:MM")
                    message = "($timeNow)Converged for seed= $seed, V=$V, α=$alpha, U=$U in $count iterations in $endTime s ($(round(endTime/60, digits = 2)) mins)\n"
                    write_to_file(message, logFileName)
                    write_to_file("----------------------------------------------------------------------------------\n", logFileName) 
                else 
                    timeNow = Dates.format(now(), "HH:MM")
                    message = "($timeNow)WARNING : (U=$U) Reached limit of Max Iteration= $maxCount\n"
                    write_to_file(message, logFileName)
                end 
                avg_delta = sum(deltaFinal)/nSites
                deltaDictV[V] = avg_delta
            end # V for loop
            deltaDictα[alpha] = deltaDictV
        end # alpha for loop
        deltaDictSeed[seed] = deltaDictα
    end # seed for loop
    newFileName = string(saveInFolder,"OP_U","_", U,".", fileFormat)
    save_file(deltaDictSeed, newFileName)
end  # main end

timestart = time()
if nnHam.type == "prebuilt"
    @everywhere begin
        vDictName = string(dataSetFolder,"vUncorrelatedDict$(N)", ".", fileFormat)
        vDict = load_file(vDictName)
        tMat = readdlm(nnHam.fileName)
    end
elseif nnHam.type == "!prebuilt"
    @everywhere begin
        _, _, nnMap , df = create_nnmat_df(N, lattice, bc) #nnMap[site][2] = nearest neighbors
        if uncorrelated == true
            vDictName = string(dataSetFolder,"vUncorrelatedDict$(N)", ".", fileFormat)
        elseif correlated == true
            vDictName = string(dataSetFolder,"vDict_$lattice$N", ".", fileFormat)
        end # if block end
        vDict = load_file(vDictName)
        tMat = create_t_matrix(N, lattice, bc, nnMap, saveAsText=tMatSave, fileName=tMatFileName)
    end
else
    throw("ERROR : choose nnHam.type to be prebuilt or !prebuilt")
end # if block end

pmap(main, UList)
elapsed = round(time() - timestart, digits = 2)
timeNow = Dates.format(now(), "HH:MM")
println("($timeNow)The elapsed time : $elapsed secs ($(round(elapsed/60, digits = 2)) mins)")
