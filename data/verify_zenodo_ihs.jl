# julia --project=. data/verify_zenodo_ihs.jl /path/to/covid19inputoutput.zip
#
# Rebuild the bundled IHS classification from Zenodo v3's processed
# `A.essential3.csv` and verify every value and bundled label. This is manual:
# the pinned source archive is 560 MB and is intentionally not a test fixture.
using DelimitedFiles
using ESRIcascade
using SHA

const _ZENODO_ENTRY = "covid19inputoutput/data/IHS_matrices_processed/A.essential3.csv"
const _ZENODO_ENTRY_SHA256 = "fe30f4ccc68493fe8e773846ca2b2134e0e4d11d7355b88c9990bd936c4aab22"
const _BUNDLED_LABELS_SHA256 = "29cbacaf97a3cb51c9e0b92603a024f3b1f9f6a45392c0c675a4bc07142ac61b"
const _BUNDLED_CLASSIFICATION_SHA256 = "5f90f03db4af0c144e989bb46d21db016ef7841bf5d9ff4e1c23de79a97f3a27"

const _IHS_DIVISION_GROUPS = Dict(
    lpad(string(division), 2, '0') => group
    for (group, divisions) in (
        "A01" => 1:1, "A02" => 2:2, "A03" => 3:3, "B" => 5:9,
        "C10-C12" => 10:12, "C13-C15" => 13:15, "C16" => 16:16,
        "C17" => 17:17, "C18" => 18:18, "C19" => 19:19, "C20" => 20:20,
        "C21" => 21:21, "C22" => 22:22, "C23" => 23:23, "C24" => 24:24,
        "C25" => 25:25, "C26" => 26:26, "C27" => 27:27, "C28" => 28:28,
        "C29" => 29:29, "C30" => 30:30, "C31_C32" => 31:32, "C33" => 33:33,
        "D35" => 35:35, "E36" => 36:36, "E37-E39" => 37:39, "F" => 41:43,
        "G45" => 45:45, "G46" => 46:46, "G47" => 47:47, "H49" => 49:49,
        "H50" => 50:50, "H51" => 51:51, "H52" => 52:52, "H53" => 53:53,
        "I" => 55:56, "J58" => 58:58, "J59_J60" => 59:60, "J61" => 61:61,
        "J62_J63" => 62:63, "K64" => 64:64, "K65" => 65:65, "K66" => 66:66,
        "L68" => 68:68, "M69_M70" => 69:70, "M71" => 71:71, "M72" => 72:72,
        "M73" => 73:73, "M74_M75" => 74:75, "N" => 77:82, "O84" => 84:84,
        "P85" => 85:85, "Q" => 86:88, "R_S" => 90:96, "" => 97:99,
    ) for division in divisions
)

function ihs_group(code::String)
    code == "9999" && return nothing
    return get(_IHS_DIVISION_GROUPS, code[1:2]) do
        error("no IHS group for NACE code $code")
    end
end

archive = only(ARGS)
isfile(archive) || error("archive not found: $archive")
source_text = read(`unzip -p $archive $_ZENODO_ENTRY`, String)
bytes2hex(sha256(source_text)) == _ZENODO_ENTRY_SHA256 || error("unexpected Zenodo entry")

raw = readdlm(IOBuffer(source_text), ',', String)
source_labels = String.(raw[1, 2:end])
source = parse.(Float64, raw[2:end, 2:end])
size(source) == (55, 55) || error("unexpected A.essential3 dimensions")
source_index = Dict(label => index for (index, label) in pairs(source_labels))

codes = ihs_industry_codes()
categories = ihs_group.(codes)
source_rows = [isnothing(category) ? 0 : source_index[category] for category in categories]
expected = fill(UInt8(1), length(codes), length(codes)) # `9999` is non-essential by construction.
for customer in eachindex(codes), supplier in eachindex(codes)
    row, column = source_rows[supplier], source_rows[customer]
    row == 0 || column == 0 || (expected[supplier, customer] = round(UInt8, 2 * source[row, column]))
end

bytes2hex(sha256(join(codes, "\n") * "\n")) == _BUNDLED_LABELS_SHA256 || error("unexpected bundled labels")
bytes2hex(sha256(vec(expected))) == _BUNDLED_CLASSIFICATION_SHA256 || error("unexpected reconstructed matrix")
expected == ihs_input_classification() || error("bundled matrix differs from Zenodo reconstruction")
println("Verified Zenodo v3 IHS reconstruction: 616 × 616 values and label order match.")
