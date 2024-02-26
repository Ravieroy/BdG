include("../src/bdg_utilities.jl")
using .BdgUtilities 

include("../src/generalutils.jl")
using .GeneralUtils

include("../scripts/params.jl")

# Create dataframes only if nnHam is not prebuilt 
if nnHam.type == "!prebuilt"
    saveDf = true
    dfName = string(dataSetFolder,"df_$lattice$N", ".csv")
    unitcell, lat, nnMap , df = create_nnmat_df(N, lattice, bc, saveAsCsv=saveDf, fileName=dfName)

else
    dfName = "NA"
end

if correlated == true
    vDictName = string(dataSetFolder,"vDict_$lattice$N", ".", fileFormat)
    vDict = make_random_V_dict(seedList, alphaList, N, df)
    vCorrDict = get_corr_V_dict(seedList, alphaList, vVals, N, df);
    save_file(vCorrDict, vDictName; key = "data")
    fileName = vDictName
elseif partial == true
    subLatArray = ["A", "B", "C"]
    newLatArray = []
    for val in 1: length(partialDisorderOn)
        if partialDisorderOn[val] == true
            push!(newLatArray, subLatArray[val])
        end 
    end 
    name = join(newLatArray)
    vPartialDictName = string(dataSetFolder,"vPartialDict_$name$N", ".", fileFormat)
    vPartialDict = create_partial_disorder_dict(partialDisorderOn, vVals, seedList, unitcell, lat, N)
    save_file(vPartialDict, vPartialDictName)
    fileName = vPartialDictName
elseif sublattice == true
    # r1 = round(subLatticeDisorder.input[2][1], digits=0)
    # r2 = round(subLatticeDisorder.input[2][2], digits=0)
    r1 = subLatticeDisorder.input[2][1]
    r2 = subLatticeDisorder.input[2][2]
    s = join(subLatticeDisorder.sites)
    vPartialDictName = string(dataSetFolder,"vSubPartialDict_$(s)_$(r1)_$(r2)_$(N)", ".", fileFormat)
    numDisorderedSitesDictName = string(dataSetFolder,"numDisorderedSitesDict_$(s)_$(r1)_$(r2)_$(N)", ".", fileFormat)
    vPartialDict, numDisorderedSitesDict = create_sublattice_disorder_dict(subLatticeDisorder, vVals, seedList, unitcell, lat, N)
    save_file(vPartialDict, vPartialDictName)
    save_file(numDisorderedSitesDict, numDisorderedSitesDictName)
    fileName = vPartialDictName

elseif uncorrelated == true    
    vUncorrelatedDict = create_uncorrelated_disorder(N, lattice, seedList, vVals, pyRandom=true, nLayer=nLayer) 
    vUncorrelatedDictName = string(dataSetFolder,"vUncorrelatedDict$(N)", ".", fileFormat)
    save_file(vUncorrelatedDict, vUncorrelatedDictName)
    fileName = vUncorrelatedDictName
else
    throw("""correlated, partial, sublattice are false.
          For uncorrelated disorder on all the sites, use
          partial=true, partialDisorderOn=[true, true, true]
          or uncorrelated=true
          """)
end 

println("""
        Summary:
        lattice = $lattice
        nLayers = $nLayer
        N = $N
        bc = $bc
        partial = $partial
        correlated = $correlated
        sublattice = $sublattice
        uncorrelated = $uncorrelated
        dataframe = $dfName
        dataset Folder = $dataSetFolder 
        filename = $fileName
        """)



