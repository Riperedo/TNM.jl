# Experimental Dataset Ingestion and Batch Loading Engine

using DelimitedFiles

"""
    load_experimental_file(filepath::AbstractString; 
                           temp_col::Int=1, 
                           cp_col::Int=2, 
                           protocol::Union{Nothing, AbstractThermalProtocol}=nothing, 
                           is_normalized::Bool=true,
                           comments::Bool=true,
                           comment_char::Char='#')

Loads experimental DSC/TMDSC measurement data from a text, `.dat`, or `.csv` file.
Ignores lines beginning with `comment_char` (such as headers `# prop1 prop2 ...`).

# Returns
- `ExperimentalDSCData` struct containing experimental temperatures, measured heat capacities, and associated protocol.
"""
function load_experimental_file(filepath::AbstractString; 
                               temp_col::Int=1, 
                               cp_col::Int=2, 
                               protocol::Union{Nothing, AbstractThermalProtocol}=nothing, 
                               is_normalized::Bool=true,
                               comment_char::Char='#')
    isfile(filepath) || throw(ArgumentError("File not found: $filepath"))

    # Read lines manually to filter comment characters
    lines = readlines(filepath)
    filtered_lines = filter(l -> !isempty(strip(l)) && !startswith(strip(l), comment_char) && !startswith(strip(l), '%'), lines)
    isempty(filtered_lines) && throw(ArgumentError("No data rows found in file: $filepath"))

    # Parse matrix using readdlm on joined lines
    raw_text = join(filtered_lines, "\n")
    data_mat = readdlm(IOBuffer(raw_text))

    size(data_mat, 2) >= max(temp_col, cp_col) || throw(DimensionMismatch("File does not contain columns $temp_col and $cp_col"))

    temps = Float64.(data_mat[:, temp_col])
    cps = Float64.(data_mat[:, cp_col])

    # Infer linear ramp protocol if none provided
    actual_proto = if protocol !== nothing
        protocol
    else
        # Default fallback: assume 20 K/min ramp across temperature range
        rate_k_s = temps[end] >= temps[1] ? (20.0 / 60.0) : (-20.0 / 60.0)
        LinearRampProtocol(temps[1], temps[end], rate_k_s)
    end

    return ExperimentalDSCData(temps, cps, actual_proto; is_normalized=is_normalized)
end

"""
    load_experimental_directory(dirpath::AbstractString; 
                                pattern::Regex=r"\\.(dat|csv|txt)\$",
                                temp_col::Int=1, 
                                cp_col::Int=2,
                                is_normalized::Bool=true)

Scans a directory for experimental files matching `pattern` and loads them into a `Vector{ExperimentalDSCData}`.
"""
function load_experimental_directory(dirpath::AbstractString; 
                                    pattern::Regex=r"\.(dat|csv|txt)$",
                                    temp_col::Int=1, 
                                    cp_col::Int=2,
                                    is_normalized::Bool=true)
    isdir(dirpath) || throw(ArgumentError("Directory not found: $dirpath"))
    files = filter(f -> occursin(pattern, f), readdir(dirpath; join=true))
    
    datasets = ExperimentalDSCData[]
    for f in files
        try
            push!(datasets, load_experimental_file(f; temp_col=temp_col, cp_col=cp_col, is_normalized=is_normalized))
        catch e
            @warn "Could not parse experimental file $f: $e"
        end
    end
    return datasets
end
