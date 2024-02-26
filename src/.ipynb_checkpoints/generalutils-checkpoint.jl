module GeneralUtils
    # Libraries required
    using JLD2, FileIO
    using HDF5, JLD
    using NPZ
    using PyCall
    using PyPlot

    @pyimport pickle

    #Module exports
    export check_rel_tol
    export save_file
    export load_file
    export write_to_file
    export show_all
    export show_matrix
    """
        check_rel_tol(oldList, newList; tol=1e-5, nRound=10)
    Returns true/false if the two lists are within tolerance value
    """
    function check_rel_tol(oldList, newList; tol = 1e-5, nRound = 10)
        tolList = tol * ones(Float64, length(newList))
        relVals = (newList - oldList) ./ newList
        relTolList = [round(abs(i), digits = nRound) for i in relVals]
        flag = relTolList <= tolList
        return flag
    end # check_rel_tol end

    """
        save_file(object, fileName; key="data")
    Saves the given object in desired format(jld, jld2, pkl, npz)
    """
    function save_file(object, fileName; key = "data")
        if last(fileName, 3) == "jld"
            save(fileName, key, object)
        elseif last(fileName, 3) == "npz"
            npzwrite(fileName, object)
        elseif last(fileName, 4) == "jld2"
            save(fileName, key, object)
        elseif last(fileName, 3) == "pkl"
            f = open(fileName, "w")
            pickle.dump(object, f, protocol = pickle.HIGHEST_PROTOCOL)
            close(f)
        else
            throw("ERROR : Possibly wrong format ~ Try jld, jld2, npz or pkl")
        end

    end # save_file end

    """
        load_file(fileName; key="data")
    Loads the file from formats(npz, jld, jld2, pkl)
    """
    function load_file(fileName; key = "data")
        if last(fileName, 3) == "npz"
            mat = npzread(fileName)
        elseif last(fileName, 3) == "jld" || last(fileName, 4) == "jld2"
            mat = load(fileName)[key]
        elseif last(fileName, 3) == "pkl"
            # load the pickle file.
            f = open(fileName, "r")
            mat = pickle.load(f)
        else
            println("ERROR : $fileName not found")
        end
    end # load_file end
    
    """
        write_to_file(message, fileName; mode="a")
    Writes the message into the file and save it locally.
    """
    function write_to_file(message, fileName; mode="a")
        f = open(fileName, mode)
        write(f, message)
        close(f)
    end


    """
        show_all(obj)
    shows the obj without any truncation
    """
    function show_all(obj)
        return show(stdout, "text/plain", obj)
    end

    function show_matrix(mat)
        PyPlot.gray()
        imshow(mat,interpolation="none")
        colorbar()
    end
end #module end
