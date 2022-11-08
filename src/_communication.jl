### Loading paths for wrapper and TS itself ###################################

const TSwrapper   = Ref("")
const TS          = Ref("")
const WrapperPy   = PyNULL()
const computeOpac = PyNULL()
const SlurmEnv    = Ref(false)

"""Load the path of the TS codes."""
load(ts_help::Ref{String}, path::String) = begin
    ts_help[] = isdir(abspath(path)) ? abspath(path) : ""
    ts_help[]
end

load_wrapper(path=joinpath(dirname(@__FILE__), "../../TurboSpectrum-Wrapper/") ) = load(TSwrapper, path)
load_TS(     path=joinpath(dirname(@__FILE__), "../../Turbospectrum_NLTE/") )    = load(TS, path)


### General Convenience #######################################################

"""The given path is converted from relative to the wrapper to absolute."""
macro inWrapper(relative_path)
	relative_path_l = esc(relative_path)
	:(_in_wrapper($(relative_path_l)))
end

"""The given path is converted from relative to the TS root dir to absolute."""
macro inTS(relative_path)
	relative_path_l = esc(relative_path)
	:(_in_ts($(relative_path_l)))
end

_in_wrapper(relative_path) = TSwrapper[] == "" ? abspath(joinpath(load_wrapper(), relative_path)) : abspath(joinpath(TSwrapper[], relative_path))
_in_ts(relative_path)      = TS[] == ""        ? abspath(joinpath(load_TS(), relative_path))      : abspath(joinpath(TS[], relative_path))

function import_wrapper(mod="source")
    if (TSwrapper[] == "") | (!(isdir(TSwrapper[])))
        error("Could not import TS wrapper from TSWrapper path. Please call load_TS() with the appropriate location first.")
    end

    copy!(WrapperPy,   _get_help_py(mod, TSwrapper[]))
    copy!(computeOpac, _get_help_py("compute_opac", TSO.@inWrapper(mod)))

    WrapperPy, computeOpac
end


### Importing Python scripts ##################################################

macro fromWrapper(mod, dir=dirname(TSO.@inWrapper("source")))
    mod_e    = esc(mod)
    path_esc = esc(dir)
    mod_s = :($mod)
    :($(mod_e) = _get_help_py($(QuoteNode(mod_s)), $path_esc))
end

macro fromTS(mod, dir=dirname(@inTS("source")))
    mod_e    = esc(mod)
    path_esc = esc(dir)
    mod_s = :($mod)
    :($(mod_e) = _get_help_py($(QuoteNode(mod_s)), $path_esc))
end


### Writing colunms in Stagger format #########################################

"""Write a column to file using the Stagger format as specified in TS"""
function write_as_stagger(name::String; teff, logg, FeH, T, rho, z=zeros(length(T)), id="TSO.jl-column")
    teff = @sprintf "%.1f" teff
    logg = @sprintf "%.5f" logg
    FeH  = @sprintf "%.2f" FeH
    n    = @sprintf "%i" length(z)

    content  = "Stagger  $(id)\n"
    content *= "$(teff) $(logg) $(FeH) $(n)\n"

    for row in eachindex(z)
        content *= @sprintf "%13.5e %13.2f %13.5e\n" z[row] T[row] rho[row]
    end

    open(name, "w") do f
        write(f, content)
    end

    nothing
end

write_as_stagger(name::String, T::Vector, rho::Vector) = write_as_stagger(name, T=T, rho=rho, z=zeros(length(T)), teff=0.0, logg=0.0, FeH=0.0)



### Communication with the TS wrapper #########################################

function prepare_TS!(setup, atmFile)
    #@info "preparing $(atmFile)"
    atmos = computeOpac.model_atmosphere()
    atmos.read(atmFile, setup.atmos_format)

    if strip(setup.atmos_format) == "marcs'"
        setup.ts_input["MARCS-FILE"] = ".true."
        atmos.path = atmFile
    elseif strip(setup.atmos_format) == "stagger"
        setup.ts_input["MARCS-FILE"] = ".false."
        atmos.path = atmFile
    else
        atmos.path = "$(setup.cwd)/$(atmos.ID).dat"
        computeOpac.write_atmos_m1d4TS(atmos,  atmos.path)
    end

    modelOpacFile = "$(setup.cwd)_$(atmos.id)_$(setup.jobID)"
    return atmos, modelOpacFile
end

function babsma!(setup, elementalAbundances; timeout="1:00:00")
    jobs = Dict()
    mem  = slurm_setup()

    # Run all job steps
    for (i, atmFile) in enumerate(setup.atmos_list)
        atmos, modelOpacFile = prepare_TS!(setup, atmFile)
        
        id  = "$(atmos.id)_$(setup.jobID)"
        job = run_babsma(setup.ts_input, elementalAbundances, atmos, modelOpacFile, id, 
                            quite=setup.debug,
                            memMB=mem, timeout=timeout)
        
        @info "running job $(id)"

        @assert !(id in keys(jobs))
        jobs[id] = job
    end

    # Wait for completion
    for id in keys(jobs)
        s = success(jobs[id])
        
        @info "   --> Babsma Job $id finished. Success: $s"
    end
end

function bsyn!(setup, elementalAbundances; timeout="1:00:00")
    jobs = Dict()
    mem  = slurm_setup()
    job_input = Dict()

    # Run all job steps
    for (i, atmFile) in enumerate(setup.atmos_list)
        atmos, modelOpacFile = prepare_TS!(setup, atmFile)

        setup.ts_input["MULTIDUMP"] = ".true.'"
        specResultFile = "./$(atmos.ID)"

        id  = "$(atmos.id)_$(setup.jobID)"
        job = run_bsyn(setup.ts_input, elementalAbundances, atmos, modelOpacFile,
                        specResultFile, id, 
                        quite=setup.debug, memMB=mem, timeout=timeout)
        
        @info "...running job $(id)"

        @assert !(id in keys(jobs))
        
        jobs[id]      = job
        job_input[id] = deepcopy([setup.ts_input, elementalAbundances, atmos, modelOpacFile, specResultFile, id])
    end

    # Wait for completion
    for id in keys(jobs)
        s = success(jobs[id])
        
        @info "   --> Bsyn Job $id finished. Success: $s"
        if !s
            @warn "One of the jobs failed. Try re-submitting..."
            job = srun_bsyn(job_input[job]..., quite=setup.debug, memMB=mem, timeout=timeout)
            s   = success(job)

            @info "   --> RE: Bsyn Job $id finished. Success: $s"
        end
        #!s && error("Job $id failed. Please check the log file in the TS folder. Exiting.")
    end
end

"""
    Creates input for the babsma.f routine and executes it
    babsma.f routine computes opacities for the give model atmosphere
    which are then used by the bsyn.f routine

    Parameters
    ----------
    ts_input : dict
        contains TS input flags
        must include the following flags:
            'MARCS-FILE'('.true.' or '.false.'),
            'ts_root' (path to TS executables bsyn.f and babsma.f)
    atmos : model_atmosphere
        for which model atmosphere to compute the opacities
    modelOpacFile : str
        where to store computed opacities
    quite : boolean
        controls details printout of the progress info
"""
function run_babsma(ts_input, elementalAbundances, atmos, modelOpacFile, id; quite=true,
                        memMB=70000/40, timeout="5:00:00")
    lmin  = @sprintf "%.3f" ts_input["LAMBDA_MIN"]
    lmax  = @sprintf "%.3f" ts_input["LAMBDA_MAX"]
    lstep = @sprintf "%.3f" ts_input["LAMBDA_STEP"]
    feh   = @sprintf "%.3f" atmos.feh

    babsma_conf = """
    'MODELINPUT:'    '$(atmos.path)'
    'LAMBDA_MIN:'    '$(lmin)'
    'LAMBDA_MAX:'    '$(lmax)'
    'LAMBDA_STEP:'   '$(lstep)'
    'MARCS-FILE:' '$(ts_input["MARCS-FILE"])'
    'MODELOPAC:' '$(modelOpacFile)'
    'METALLICITY:'    '$(feh)'
    'HELIUM     :'    '0.00'
    'R-PROCESS  :'    '0.00'
    'S-PROCESS  :'    '0.00'
    """

    babsma_conf = babsma_conf *
            "'INDIVIDUAL ABUNDANCES:'   '$(length(elementalAbundances))' \n"

    for i in eachindex(elementalAbundances)
        z, abund = elementalAbundances[i]
        
        zs     = @sprintf "%.0f" z
        abunds = @sprintf "%5.3f" abund

        babsma_conf = babsma_conf * " $(zs) $(abunds) \n"
    end

    ddir    = @inTS("")
    cdir    = pwd()
    command = if SlurmEnv[] 
        !isnothing(timeout) ?
            Cmd(`srun -N 1 -n 1 -c 1 --mem-per-cpu=$(memMB)M --time=$(timeout) --exclusive -D $(ddir) ./exec/babsma_lu`) :
            Cmd(`srun -N 1 -n 1 -c 1 --mem-per-cpu=$(memMB)M --exclusive -D $(ddir) ./exec/babsma_lu`)
    else
        cd(ddir)
        `./exec/babsma_lu`
    end

    p = Pipe()
    r = run(pipeline(command, stdout=joinpath(ddir, "babsma_$(id).log"), 
                              stderr=joinpath(ddir, "babsma_$(id).err"), 
                              stdin =p), wait=false)
    write(p, babsma_conf)
    close(p)

    !SlurmEnv[] && cd(cdir)

    return r
end

"""
    Creates input for the bsyn.f routine and executes it
    bsyn.f runs spectral synthesis based on the opacities
    computed previsously by babsma.f

    Parameters
    ----------
    ts_input : dict
        contains TS input flags
        must include the following flags:
            'NLTE' ('.true.' or '.false.'),
            'LAMBDA_MIN', 'LAMBDA_MAX', 'LAMBDA_STEP',
            'MARCS-FILE'('.true.' or '.false.'),
            'NFILES' (how many linelists provided, integer),
            'LINELIST' (separated with new line),
            'ts_root' (path to TS executables bsyn.f and babsma.f)
    elementalAbundances : list
        contains atomic numbers and abundances of the elements requested for spectral synthesis,
        e.g.
        [ [26, 7.5], [8, 8.76] ]
    atmos : model_atmosphere
        for which model atmosphere to compute the opacities
    modelOpacFile : str
        path to the file storing opacities computed by babsma
        returned by compute_babsma
    specResultFile : str
        where to save computed spectrum
    nlteInfoFile : str
        path to the configuration file that controls inclusion of NLTE to TS
        returned by create_NlteInfoFile
        if None, spectrum will be computed in LTE
    quite : boolean
        controls details printout of the progress info
"""
function run_bsyn(ts_input, elementalAbundances, atmos, modelOpacFile, specResultFile, id; quite = true,
                                                memMB=70000/40, timeout="5:00:00")
    lmin  = @sprintf "%.3f" ts_input["LAMBDA_MIN"]
    lmax  = @sprintf "%.3f" ts_input["LAMBDA_MAX"]
    lstep = @sprintf "%.3f" ts_input["LAMBDA_STEP"]

    bsyn_config = """
    'NLTE :'          '$(ts_input["NLTE"])'
    'LAMBDA_MIN:'    '$(lmin)'
    'LAMBDA_MAX:'    '$(lmax)'
    'LAMBDA_STEP:'   '$(lstep)'
    'INTENSITY/FLUX:' 'Flux'
    'MARCS-FILE:' '$(ts_input["MARCS-FILE"])'
    'MODELOPAC:'        '$(modelOpacFile)'
    'RESULTFILE :'    '$(specResultFile)'
    'HELIUM     :'    '0.00'
    'NFILES   :' '$(ts_input["NFILES"])'
    $(ts_input["LINELIST"])
    """
    
    bsyn_config = bsyn_config * " 'SPHERICAL:'  '.false.' \n"
    bsyn_config = bsyn_config * """
    30
    300.00
    15
    1.30
    """

    bsyn_config = bsyn_config *
            "'INDIVIDUAL ABUNDANCES:'   '$(length(elementalAbundances))' \n"

    for i in eachindex(elementalAbundances)
        z, abund = elementalAbundances[i]
        
        zs     = @sprintf "%.0f" z
        abunds = @sprintf "%5.3f" abund

        bsyn_config = bsyn_config * " $(zs) $(abunds) \n"
    end

    ddir    = @inTS("")
    cdir    = pwd()
    command = if SlurmEnv[] 
        !isnothing(timeout) ?
            `srun -N 1 -n 1 -c 1 --mem-per-cpu=$(memMB)M --time=$(timeout) --exclusive -D $(ddir) ./exec/bsyn_lu` :
            `srun -N 1 -n 1 -c 1 --mem-per-cpu=$(memMB)M --exclusive -D $(ddir) ./exec/bsyn_lu`
    else
        cd(ddir)
        `./exec/bsyn_lu`
    end
    
    p = Pipe()
    r = run(pipeline(command, stdout=joinpath(ddir, "bsyn_$(id).log"), 
                              stderr=joinpath(ddir, "bsyn_$(id).err"), 
                              stdin =p), wait=false)
    write(p, bsyn_config)
    close(p)

    !SlurmEnv[] && cd(cdir)

    return r
end



### Output handling ###########################################################

function move_output(move_to=TSO.@inWrapper("example/garbage"), wrapper_path=pwd(), ts_path=TSO.@inTS(""))
    if !isdir(move_to)
        mkdir(move_to)
    else
        rm(move_to, recursive=true)
        mkdir(move_to)
    end
    @assert isdir(move_to)

    # Clean the models from the wrapper
    @info "Moving files from $(wrapper_path) to $(move_to)"

    paths = glob("_TSOeos*", wrapper_path)
    for path in paths
        goal = joinpath(move_to, basename(path))
        mv(path, goal)
    end

    paths = glob("TSOeos_*.dat", wrapper_path)
    for path in paths
        goal = joinpath(move_to, basename(path))
        mv(path, goal)
    end

    # Clean the models from the wrapper
    @info "Moving files from $(ts_path) to $(move_to)"

    paths = glob("TSOeos*.multi", ts_path)
    for path in paths
        goal = joinpath(move_to, basename(path))
        mv(path, goal)
    end

    # Babsma logs
    paths = glob("babsma*.log", ts_path)
    for path in paths
        goal = joinpath(move_to, basename(path))
        mv(path, goal)
    end
    paths = glob("babsma*.err", ts_path)
    for path in paths
        goal = joinpath(move_to, basename(path))
        mv(path, goal)
    end

    # Bsyn logs
    paths = glob("bsyn*.log", ts_path)
    for path in paths
        goal = joinpath(move_to, basename(path))
        mv(path, goal)
    end
    paths = glob("bsyn*.err", ts_path)
    for path in paths
        goal = joinpath(move_to, basename(path))
        mv(path, goal)
    end
end



### Read slurm setup ##########################################################

slurm_setup() = begin
    SlurmEnv[] = if "SLURM_NTASKS" in keys(ENV)
        true
    else
        @warn "No slurm environment detected. Jobs will be submitted without resource management!"
        false
    end

    mem = if !("SLURM_NTASKS_PER_NODE" in keys(ENV))
        1000
    else
        @assert "SLURM_NTASKS_PER_NODE" in keys(ENV)
        @assert "SLURM_MEM_PER_NODE" in keys(ENV)

        tasks, m = parse.(Int, [ENV["SLURM_NTASKS_PER_NODE"], ENV["SLURM_MEM_PER_NODE"]])

        @info "Slurm setup: \n  tasks per note=$(tasks)\n  mem per task (MB)=$(floor(Int, m / tasks))"

        floor(Int, m / tasks)
    end

    mem
end

###############################################################################