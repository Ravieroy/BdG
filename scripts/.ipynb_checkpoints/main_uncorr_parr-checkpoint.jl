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

@everywhere function main(seed)
    #------Logging part-----------
    logFileName = string(logFolder, "info","_seed_$seed", "_id_$(myid())", ".log")
    dateFormat = "yyyy-mm-dd HH:MM:SS"
    startDate = Dates.format(now(), dateFormat)
    message1 = "\nBdG Calculation for N=$N for seed=$seed started : $startDate\n"
    message2 = string("number of processes : $(nprocs())\n", "number of workers : $(nworkers())\n" )
    message = string(message1, message2, "\n")
    write_to_file(message, logFileName)
    #-----------Logging part end ------

    deltaDict = Dict()
    nAvgDict = Dict()
    eGapDict = Dict()
    deltaGapDict = Dict()
    evecsDict = Dict()
    evalsDict = Dict()
    for alpha in alphaList
        deltaDictV = Dict()
        nAvgDictV = Dict()
        eGapDictV = Dict()
        deltaGapDictV = Dict()
        evecsDictV = Dict()
        evalsDictV = Dict()
        for V in vVals
            #------Logs-------
            timeNow = Dates.format(now(), "HH:MM")
            message = "Running($timeNow) : seed=$seed V=$V \n"
            write_to_file(message, logFileName)
 

            vList = vDict[seed][V]
            deltaFinal, nAvgFinal, eGap, deltaGap, evecs, evals, count, endTime, isConverged = 
                 run_self_consistency_numpy(deltaOld, mu, N, nAvgOld, nExp, lattice, seed, tMat, U, vList, tol=tol, nLayer=nLayer)
            if isConverged == true 
                message = "Converged for seed= $seed, V=$V, alpha=$alpha in $count iterations in $endTime s ($(round(endTime/60, digits = 2)) mins)\n"
                write_to_file(message, logFileName)
                write_to_file("---------------------------------------------------------------------------\n", logFileName)
            end
            deltaDictV[V] = deltaFinal
            nAvgDictV[V] = nAvgFinal
            eGapDictV[V] = eGap
            deltaGapDictV[V] = deltaGap
            evecsDictV[V] = evecs
            evalsDictV[V] = evals
        end #V loop end
        deltaDict[alpha] = deltaDictV
        nAvgDict[alpha] = nAvgDictV
        eGapDict[alpha] = eGapDictV
        deltaGapDict[alpha] = deltaGapDictV
        evecsDict[alpha] = evecsDictV
        evalsDict[alpha] = evalsDictV
    end # alpha loop end  
     # code block to save dictionary locally
    dictList = [deltaDict, nAvgDict, eGapDict, deltaGapDict, evecsDict, evalsDict]
    for (key, value) in store
        if value[1] == true
            baseFileName = value[2]
            newFileName = string(saveInFolder,"/", baseFileName,"_", seed,".", fileFormat)
            save_file(dictList[value[3]], newFileName)
        end 
    end # store loop end
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

        vUncorrelatedDictName = string(dataSetFolder,"vUncorrelatedDict$(N)", ".", fileFormat) 
        vDict = load_file(vUncorrelatedDictName)
        tMat = create_t_matrix(N, lattice, bc, nnMap, saveAsText=tMatSave, fileName=tMatFileName)
    end  
else
    throw("ERROR : choose nnHam.type to be prebuilt or !prebuilt")
end # if block end 

pmap(main, seedList)
elapsed = round(time() - timestart, digits = 2)
timeNow = Dates.format(now(), "HH:MM")
println("($timeNow)The elapsed time : $elapsed secs ($(round(elapsed/60, digits = 2)) mins)") 
