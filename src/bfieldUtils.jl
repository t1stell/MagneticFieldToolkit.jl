
"""
vector2range(v)

Convert an evenly spaced vector points to a range
"""
function vector2range(v::AbstractVector)
  if typeof(v) <: AbstractRange
    return v
  else
    v_diff = @view(v[2:end]) .- @view(v[1:end-1])
    v_diff_mean = sum(v_diff)/length(v_diff)
    @assert all(v_diff .≈ v_diff_mean) "Step between vector elements must be consistent"
    return range(v[1],step=v_diff_mean,length=length(v))
  end
end

function irft(f::Vector{ComplexF64}, M::Integer, x::Float64)

    res = f[1]
    N = length(f)
    for i in 2:(N - iseven(M))
       res += f[i]*exp(im*2π*(i-1)*x/M) + conj(f[i])*exp(-im*2π*(i-1)*x/M)
    end
    return real((res + iseven(M)*f[end])/M)

end
"""
irfftAtPoint(rfft, np, x)

Calculate an inverse rfft of original length np at point x
"""
function irfftAtPoint(rfft::Vector{ComplexF64}, np::Integer, x::Float64)
  ans = sum([rfft[i] * exp(im * x * (i-1) / np) for i in 1:length(rfft)])
  ms = div(np,2)
  ans += sum([conj(rfft[np-i]) * exp(im * x * (i+1) / np ) for i in ms:np-2])
  return real(ans)/np
end
