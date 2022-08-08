

struct CoilFilament{FT}
  x::Array{FT, 1}
  y::Array{FT, 1}
  z::Array{FT, 1}
  current::FT
end

struct CoilFamily{FT}
  coil::Array{CoilFilament{FT}, 1}
  name::String
end
"""
  The coil hierarchy uses the convention common in various coil formats, like
  the ones used by vmec and mgrid. 
  Each coil set is broken into a family of coil filaments that all have
  the same currents, and each family is broken up into individual filaments
  Accessing a specific coil is done by:
  coil_set.family[i].coil[j]. etc
"""
struct CoilSet{FT}
  family::Array{CoilFamily{FT},1}
end

"""
  read_vmec_coils(filename)

Reads a coil file in the vmec format and outputs a coil_set 
Right now the code assumes all coils exist and no explicit rotations
are calculated. Check coil file to verify
"""
function read_vmec_coils(filename::String)
  #we read through the file twice, first to determine how many coil families
  #there are and how many coils in each family, then again
  #to save the actual coil data
  family_names = Vector{String}(undef, 0)
  coils_per_family = Vector{Int64}(undef, 0)
  lines = readlines(filename)
  coilends = Vector{Int64}(undef, 0)
  family_labels = Vector{Int64}(undef, 0)
  for (i, ln) in enumerate(lines)
    dum = split(ln)
    if length(dum) < 5
      continue
    end
    family_name = dum[6]
    family_label = parse(Int, dum[5])
    while family_label > length(coils_per_family)
      push!(coils_per_family, 0)
      push!(family_names, "")
    end
    coils_per_family[family_label] += 1
    if family_names[family_label] == ""
      family_names[family_label] = family_name
    end
    push!(coilends, i)
    push!(family_labels, family_label) 
  end
  nfamilies = length(family_names)
  nfilaments = length(coilends)
  
  #Now read back through the file and assign stuff
  coil_set = Vector{CoilFamily{Float64}}(undef, nfamilies)
  for fam_index in 1:nfamilies
    coil_family = Vector{CoilFilament{Float64}}(undef, coils_per_family[fam_index])
    #get all the indices
    family_indices = findall(x->x==fam_index, family_labels)
    for (j, coil_index) in enumerate(family_indices)
      #guess at start index
      coil_index == 1 ? start_index = 1 : start_index = coilends[coil_index-1]+1
      end_index = coilends[coil_index]
      #we can't preassign the vectors bc there might be junk values
      xc = Vector{Float64}(undef, 0)
      yc = Vector{Float64}(undef, 0)
      zc = Vector{Float64}(undef, 0)
      current = nothing
      for j in start_index:end_index
        dum = split(lines[j])
        if length(dum) < 4
          continue
        end
        push!(xc, parse(Float64, dum[1]))
        push!(yc, parse(Float64, dum[2]))
        push!(zc, parse(Float64, dum[3]))
        if current == nothing 
          current = parse(Float64, dum[4])
        end
      end
      coil_family[j] = CoilFilament(xc, yc, zc, current)
    end
    coil_set[fam_index] = CoilFamily(coil_family, family_names[fam_index])
  end
  return CoilSet(coil_set)
  
end  
