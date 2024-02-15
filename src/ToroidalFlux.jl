function poincare_plot_at_half_radius(mgrid::Union{MagneticField{T}, CoilSet{T}},
                                      start_radius::Real;
                                      nfp=4,
                                      ϕ_step = π/100) where {T}
    mgnfp = nfp
    try
        mgnfp = mgrid.nfp
    catch
        #no nfp set
        mgnfp = nfp
    end
    nfp = mgnfp
    ϕ_half = π/nfp

    a = follow_field(mgrid, [start_radius, ϕ_half, 0.0], 1001*ϕ_half, ϕ_step = ϕ_step, poincare=true, 
                     poincare_res = 2*π/nfp)
    r = [a.u[i][1] for i in 1:length(a.u)]
    z = [a.u[i][2] for i in 1:length(a.u)]

    return r, z
end

"""
    function iota_from_poincare(r, z, nfp)

    Given a set of points r, z (AbstractArrays) this function will estimate the rotational transform, iota
    The points must not be resort from the field line following

    It does this by comparing individual points and counting the amount of times you pass across the θ=π line
    in the plane.  This will fail if iota is greater than nfp/2, probably.
"""
function iota_from_poincare(r::AbstractArray{T}, z::AbstractArray{T}, nfp::Integer) where {T}
    r0 = sum(r)/length(r)
    z0 = sum(z)/length(z)
    times_around = 0
    prev_angle = atan(z[1]-z0, r[1]-r0)

    for i in 2:length(r)
        current_angle = atan(z[i]-z0, r[i]-r0)
        if prev_angle > 0 && current_angle < 0
            times_around += 1
        end
        prev_angle = current_angle
    end
    #assume we went retrograde
    if times_around > length(r)/2
        times_around = length(r) - times_around
    end
    return times_around * nfp / length(r)
end

"""
    function iota_at_r(mgrid, radius)

    Calculate iota at a given major radius value.  See "iota_from_poincare" for how it works

"""
function iota_at_r(mgrid::Union{MagneticField{T}, CoilSet{T}},
                   radius::Real;
                   ) where {T}

    r,z = poincare_plot_at_half_radius(mgrid, radius, nfp=mgrid.nfp)
    return iota_from_poincare(r, z, mgrid.nfp)
end

function sort_poincare(r::AbstractArray{T}, z::AbstractArray{T}) where {T}
    r0 = sum(r)/length(r)
    z0 = sum(z)/length(z)
    angles = atan.(z.-z0, r.-r0)
    #make from 0 to 2π
    angles = mod.(angles, 2π)
    sort_indices = sortperm(angles)
    return r[sort_indices], z[sort_indices], angles[sort_indices]
end

function smooth_poincare(r::AbstractArray{T}, z::AbstractArray{T}; smooth::Integer=5) where {T}
    #assumes already sorted
    
    npoints = length(r)
    #append values to both ends
    rpad = vcat(r[end-smooth+1:end], r, r[1:smooth])
    zpad = vcat(z[end-smooth+1:end], z, z[1:smooth])
    rsmooth = [sum(rpad[i-smooth:i+smooth])/(2*smooth+1) for i in smooth+1:npoints+smooth]
    zsmooth = [sum(zpad[i-smooth:i+smooth])/(2*smooth+1) for i in smooth+1:npoints+smooth]
    return rsmooth, zsmooth
end

function make_interpolation(r::AbstractArray{T}, z::AbstractArray{T}) where {T}
    r,z = sort_poincare(r, z)
    r,z = smooth_poincare(r, z)
    #resort to fix problems with smoothing crossing the boundary
    r,z = sort_poincare(r, z)
    upper_half_indices = z .>= 0
    z = z[upper_half_indices]
    r = r[upper_half_indices]
    z[1] = 0.0
    z[end] = 0.0

    r = reverse(r)
    z = reverse(z)
    #have to handle the case where it's pinches in slightly on the inboard side.  Just kill these values
    mindex = findmin(r)[2]
    r = r[mindex:end]
    z = z[mindex:end]
    z[1] = 0.0

    return linear_interpolation(r,z)
end   

function toroidal_flux_integration(mgrid::Union{MagneticField{T}, CoilSet{T}}, 
                              itp::Interpolations.Extrapolation; 
                              ) where {T}
    ϕ = π/mgrid.nfp
    rmin = itp.itp.knots[1][1]
    rmax = itp.itp.knots[1][end]
    
    function flux_integrand(v::SVector{R}) where {R}
        r = v[1]
        if r < rmin || r > rmax
            return 0.0
        end
        z = v[2]*itp(r)
        return mgrid(r,ϕ,z)[2]*itp(r)
    end

    return 2*hcubature(flux_integrand, [rmin, 0.0], [rmax, 1.0], rtol=1.0E-6)[1] 
end

"""
    function toroidal_flux(mgrid, start_radius)

    The mgrid can either be a MagneticGrid object in the julia format (this can be
    loaded directly from an mgrid.nc file or a bmw output file)
    Coil sets should also work, although the calculation may be prohibitively slow

    It is expected that you are using the ϕ=π/nfp plane, or the half-period plane
    as this is up down symmetric and tends to be better behaved than the ϕ=0 plane

    Also needed is a radius to generate the Poincare plot from, this probably needs to be
    determined beforehand via some plotting

    The code does the following steps
    1: calculates a poincare plot at the desired radius
    2: sorts the points in angle order (needs to be a star domain from the geometric center)
    3: smooths the data and resorts
    4: fits an interpolation spline (this may cut off indentations, which hopefully are small errors)
    5: integrates B_ϕ over the half domain

"""
function toroidal_flux(mgrid::Union{MagneticField{T}, CoilSet{T}},
                       start_radius::Real) where {T}
    r,z = poincare_plot_at_half_radius(mgrid, start_radius)
    itp = make_interpolation(r, z)
    return toroidal_flux_integration(mgrid, itp)
end
