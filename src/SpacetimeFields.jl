module SpacetimeFields

using NetCDF, DataFrames, Statistics
ncm = NetCDF
import Base.convert, Base.copy, Base.show

export extent, stfield
export copy, llgrid, nc2field, convert, show
export array2field

"""
    extent(xmin, xmax, ymin, ymax)

An `extent` object
"""
struct extent
    xmin::Float64
    xmax::Float64
    ymin::Float64
    ymax::Float64
end


"""
    stfield(data, lon, lat, good, time)

A `stfield` object
"""
mutable struct stfield
    data::Array
    lon::Vector
    lat::Vector
    good::BitArray
    time::Vector
end


"""
    copy(x::stfield; layers::AbstractVector=[1])

Creates a copy of an `stfield` object
"""
function copy(x::stfield; layers::AbstractVector=[1])
    m1 = x.data[:,:,layers]
    stfield(m1, x.lon, x.lat, x.good, collect(layers))
end


"""
    llgrid(fname::String; start=[1], count=[-1])

Read a NetCDF file with the given `fname`, and return  a grid of longitude and latitude values 
"""
function llgrid(fname::String; start=[1], count=[-1] )
    lon = ncread(fname,"lon", start, count)
    lat = ncread(fname, "lat", start, count)
    ll = vcat([ [i j ] for i in lon, j in lat ]...)
    lli = vcat([ [i j ] for i in 1:length(lon), j in 1:length(lat) ]...)
    return ll, lli
end


"""
    llgrid(x::stfield)

Return a grid of longitude and latitude values from a `stfield` object
"""
function llgrid(x::stfield)
    ll = vcat([ [i j ] for i in x.lon, j in x.lat ]...)
    lli = vcat([ [i j ] for i in 1:length(x.lon), j in 1:length(x.lat) ]...)
    return ll, lli
end


"""
    nc2field(fname::String, varname::String, ext::extent; mean_dims::AbstractVector)

Convert NetCDF data to a `stfield` object, after averaging over dimensions `mean_dims` (optional)
"""
function nc2field(fname::String, varname::String, ext::extent; mean_dims::AbstractVector=[])
    ll,lli = llgrid(fname)
    imin = argmin( hypot.(ll[:,1].-ext.xmin, ll[:,2].-ext.ymin ) )
    imax = argmin( hypot.(ll[:,1].-ext.xmax, ll[:,2].-ext.ymax ) )

    a = lli[imin,:]
    b = lli[imax,:]
    d = abs.(b - a) .+ [1,1]
    sign(b[2]-a[2])==-1 && (a[2] = b[2])
    ncf = ncm.open(fname)
    start = ncm.defaultstart(ncf[varname])
    count = -ncm.defaultstart(ncf[varname])
    start[1:2] = a; count[1:2] = d
    a1 = ncm.readvar(ncf[varname], start=start, count=count)
    !isempty(mean_dims) && (a1 = mean(a1, dims=mean_dims))
    a1 = dropdims(a1, dims=(findall(size(a1).== 1)...,))
    
    miss_value = ncgetatt(fname, varname,"missing_value")
    good = vec(a1[:,:,1] .≠ miss_value)
    lon = ncm.readvar(ncf["lon"], start=a[1:1], count=d[1:1])
    lat = ncm.readvar(ncf["lat"], start=a[2:2], count=d[2:2])
    time = collect(1:size(a1, 3))
    ncm.close(ncf)
    return stfield(a1, lon, lat, good, time)
end

"""
    array2field(field::AbstractArray, lon, lat, ext::extent)

Convert a 3d-array to a `stfield` object 
"""
function array2field(field::AbstractArray, lon::T, lat::T,  ext::extent) where {T<:AbstractVector}

    ll = vcat([ [i j] for i in lon, j in lat ]...)
    lli = vcat([ [i j ] for i in 1:length(lon), j in 1:length(lat) ]...)

    imin = argmin( hypot.(ll[:,1].-ext.xmin, ll[:,2].-ext.ymin ) )
    imax = argmin( hypot.(ll[:,1].-ext.xmax, ll[:,2].-ext.ymax ) )

    a = lli[imin,:];    b = lli[imax,:]
#    dump(a); dump(b)    
    
    data = field[a[1]:b[1], a[2]:b[2], :]
    lon = lon[a[1]:b[1]]
    lat = lat[a[2]:b[2]]
    good = trues(length(data[:,:,1]))
    time = collect(1:size(field, 3))

    return stfield(data, lon, lat, good, time)    
end


"""
    convert(::Type{Matrix}, x::stfield)

Convert `stfield` object to a `Matrix`.
"""
function convert(::Type{Matrix}, x::stfield)
    a, b, n = size(x.data)
    d0 = reshape(x.data, a*b, n)'
    d1 = d0[:,x.good]
    return(d1)
end


"""
    convert(::Type{DataFrame}, x::stfield, time::AbstractVector)

Convert `stfield` object to a `DataFrame`
"""
function convert(::Type{DataFrame}, x::stfield, time::AbstractVector)
    ll = vcat([[i j ] for i in x.lon, j in x.lat ]...)[x.good,:]     
    D = vcat( [DataFrame(lon=ll[:,1],lat=ll[:,2], z = x.data[:,:,k][x.good], g="$k") for k in time]...) 
    return D
end


"""
    convert(::Type{DataFrame}, x::stfield, time, label::Vector{String})

Convert `stfield` object to a `DataFrame`
"""
function convert(::Type{DataFrame}, x::stfield, time, label::Vector{String})
    ll = vcat([ [i j] for i in x.lon, j in x.lat ]...)[x.good,:]     
    D = vcat( [DataFrame(lon=ll[:,1], lat=ll[:,2], z=x.data[:,:,k][x.good], g=label[k]) for k in time]...) 
    return D
end


"""
    show(io::IO,x::stfield)

Display info about a `stfield` object
"""
function show(io::IO, x::stfield)
    println(typeof(x))
    print("Data: "); print(typeof(x.data)); println(size(x.data))
    print("Longitude: "); println(x.lon')
    print("Latitude: "); println(x.lat')
    print("Time: "); println(x.time')
end


end # module
