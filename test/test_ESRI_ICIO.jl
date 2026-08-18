include("../src/ESRIcascade.jl")
using .ESRIcascade
using DelimitedFiles

GetNace2(x::AbstractVector{<:AbstractString}) = map(x->x[1:2], x)

codesjl = ihs_industry_codes()
ICIO_P = vec(readdlm("data/ICIO/p_icio.csv",',',String))
ICIO_W = readdlm("./data/ICIO/W_ICIO.csv",',',Float64,header=false)

EMS = ihs_input_classification()
ids = map(x->findfirst(isequal(x),GetNace2(codesjl)),GetNace2(ICIO_P))

blabba = IndustryInfo(ids, EMS)
Economy = ESRIEconomy(ICIO_W, blabba)

ESRI = esri(Economy)
ESRI_R = readdlm("./data/ICIO/ESRI_icio.csv",',',Float64,header=false)
@assert all(isapprox.(ESRI,ESRI_R[:,1],atol=1e-6)) "Error: ESRI results do not match reference data"