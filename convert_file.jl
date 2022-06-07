# Run on julia 1.7.2
# First time, run `using Pkg; Pkg.activate(); Pkg.instantiate()`

using CSV
using DataFrames
using Statistics

snps = CSV.read("data/resonance_snps.csv", DataFrame)
subset!(snps, :mean_cog => ByRow(!isnan))
modeldata = CSV.read("GenotypeData.raw", DataFrame; delim=' ')

prim_allele = Dict(id => first(allele) for (id, allele) in split.(names(modeldata, r"^rs"), '_'))

for id in names(snps, r"^rs")
    haskey(prim_allele, id) || continue
    snps[!, string(id, "_", prim_allele[id])] = map(eachrow(snps)) do row
        count(a-> a == prim_allele[id], collect(row[id]))
    end
end

#-

rename!(snps, "ECHOProtocolDChild"=> "IID")
snps.FID = map(snps.IID) do id
    join(split(id, "-")[1:2], '-')
end

snps.SEX = map(snps.childGender) do g
    ismissing(g) && return 0
    g == "Male" && return 1
    g == "Female" && return 2
    throw(ArgumentError("Incorrect format: $g"))
end

snps.PAT .= 0
snps.MAT .= 0

med_cog = median(snps.mean_cog)

snps.PHENOTYPE = map(snps.mean_cog) do c
    c < med_cog ? 1 : 2
end

CSV.write("resonance_genotype.raw", select(snps, Cols(
        "FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE", r"rs\d+_[ATGC]")
    ); delim=' '
)