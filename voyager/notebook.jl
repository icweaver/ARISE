### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ a7566d2f-b355-40aa-ada7-7a8830e372b3
begin
# Analysis
using Chain, FFTW, Peaks, Statistics

# Data viz
using PlutoUI, PlutoPlotly
end

# ╔═╡ ca1c0637-52c7-4168-9138-f5b91f9578dc
md"""
# 🛰️ Voyager Lab
"""

# ╔═╡ 525b7253-59d6-4253-a61f-75db1b61dd34
md"""
## Step 1: Load data 📦

Load the data
"""

# ╔═╡ 2f45818d-34bc-4bd2-a4ac-c9552a7b144e
# str = let
# 	url = "https://www.dropbox.com/scl/fi/6s0p463op00f0vpykpvfj/Voyager_output.ci8?rlkey=t7l2c32w91h46n9p4rj3goife&e=1&st=re8jzzpu&dl=1"
# 	io = IOBuffer()
# 	download(url, io)
# 	String(take!(io))
# end

# ╔═╡ 158c0826-8eb4-4b0e-9d2b-65b6e02ebb44
filtered_data = let
	f = joinpath("data", "Voyager_output.ci8") # OS-agnostic file paths
	N = filesize(f) ÷ 8 # Int8 chunks
	read!(f, Vector{ComplexF32}(undef, N)) # Equivalent of np.complex64
end

# ╔═╡ aeef9d9a-a111-42e6-94af-b47cc9c0c667
md"""
## Step 2: Analyze frequencies 💻

Frequency Axis Computation
"""

# ╔═╡ b2c12094-bd4f-4fd8-9c79-2f2a351605c1
begin
	# Acquisition parameters
	const center_freq = 8431e6    # Center frequency in Hz
	const fs = 2.92969e6          # Sampling frequency in Hz
	const FFT_size = 16384        # Number of samples per FFT segment
end;

# ╔═╡ 282076e3-d070-438b-b30e-8474fa97a3b9
freqs_corrected = fftshift(fftfreq(FFT_size, fs)) .+ center_freq

# ╔═╡ 762d18d8-28fa-445c-b615-87189d7e5a85
md"""
**Frequency range:**
$(round(first(freqs_corrected)/1e6; digits=2)) Hz to
$(round(last(freqs_corrected)/1e6; digits=2))
"""

# ╔═╡ b2aa6bf1-9943-47db-bc24-63ba3f90446e
md"""
## Step 3: Plot 📈

Plot the Entire Power Spectrum
"""

# ╔═╡ 9b8715fb-8921-4e9f-8f46-7d68a5a69e98
# Segment the data into FFT-sized chunks
num_segments = length(filtered_data) ÷ FFT_size

# ╔═╡ f1916076-f376-420d-8a91-00039badec1c
trimmed_data = filtered_data[begin : num_segments*FFT_size] # Trim extra samples

# ╔═╡ 6c06e8a7-9565-4ed1-8369-ccbeb240e820
data_segments = reshape(trimmed_data, FFT_size, num_segments)

# ╔═╡ be7a8621-a11b-4b9b-b87f-ce15ea1d9236
power_spectra = abs.(fftshift(fft(data_segments, 1))).^2

# ╔═╡ 369e9413-d489-4a23-826c-9b14387924dc
avg_power_spectrum = 10 * log10.(mean(power_spectra, dims=2) .+ 1e-12) |> vec

# ╔═╡ edb8cf40-3105-4ce6-b312-37461fb1cd23
scatter(x=freqs_corrected/1e6, y=avg_power_spectrum) |> plot

# ╔═╡ 43e78d19-e465-43c3-9b7a-7e7e8c8c382b
md"""
## Step 4: ROI 🔎
"""

# ╔═╡ 0486c276-88fd-4dce-b991-e12ec75b699a
md"""
**Lower bound:** $(@bind roi_min NumberField(8429:0.1:8431; default=8430.1)) MHz
**Upper bound:** $(@bind roi_max NumberField(8431:0.1:8433; default=8431.9)) MHz
"""

# ╔═╡ b41854aa-29ea-4806-a12e-a3585a446539
freq_mask = @. roi_min*1e6 ≤ freqs_corrected ≤ roi_max*1e6

# ╔═╡ 8bd511ab-d7af-4f68-9a29-9aaed83b3774
freqs_zoomed = @views freqs_corrected[freq_mask]

# ╔═╡ cdbf8282-69b4-4a9d-9f58-e7aa270bf868
avg_power_spectrum_zoomed = @views 10 * log10.(mean(power_spectra[freq_mask, :], dims=2) .+ 1e-12) |> vec

# ╔═╡ 01ae457d-0093-4e92-bab4-17422ce3b97a
scatter(x=freqs_zoomed/1e6, y=avg_power_spectrum_zoomed) |> plot

# ╔═╡ 78116562-b23d-4242-b4b9-f3bb669166f0
md"""
## Step 5: Peak detection 🗻

Absolute peak
"""

# ╔═╡ 3f4eaf24-3855-4f24-8024-b17a652aa1ad
md"""
param | value | description
:--   | :--   | :--
peak\_height\_offset | $(@bind peak_height_offset NumberField(1:0.1:100; default=30)) | dB offset below the ROI's maximum power to set threshold
min_distance | $(@bind min_distance NumberField(1:100; default=10)) | Minimum separation between peaks (in index units)
prominence | $(@bind prominence NumberField(1:0.1:100; default=5)) | Minimum prominence (in dB) for a peak
"""

# ╔═╡ 960903a9-2de8-4dc7-8883-da968562889e
peak_height = maximum(avg_power_spectrum_zoomed) - peak_height_offset

# ╔═╡ 7e6b20a7-2a8d-48c4-8422-d371b81da6a5
pks = @chain avg_power_spectrum_zoomed begin
	findmaxima(min_distance)
	peakproms(min=prominence)
	peakheights(min=peak_height)
end

# ╔═╡ 78c6aebd-25b6-4890-bc6c-e8e4f63eef67
let
	p1 = scatter(x=freqs_zoomed/1e6, y=avg_power_spectrum_zoomed)

	mask = pks.indices
	p2 = scatter(x=freqs_zoomed[mask]/1e6, y=avg_power_spectrum_zoomed[mask];
		mode = :markers,
		marker = attr(symbol=:x),
	)
	
	plot([p1, p2])
end

# ╔═╡ c93a0c97-f322-41ae-9f6d-ce5647b46b2d
md"""
## Step 6: SNR computation 📢
"""

# ╔═╡ f5343371-806b-44d6-aa0d-f83444217ecf
md"""
**Quiet region 1:** $(@bind quiet_1 RangeSlider(8430:0.1:8431; default=8430.4:0.1:8430.9)) MHz

**Quiet region 2:** $(@bind quiet_2 RangeSlider(8431:0.1:8432; default=8431.2:0.1:8431.6)) MHz
"""

# ╔═╡ af6ddaff-6be3-44fa-bfc6-03e7ea1b4c32
let
	p1 = scatter(x=freqs_zoomed/1e6, y=avg_power_spectrum_zoomed)

	fig = plot(p1)

	add_vrect!(fig, extrema(quiet_1)...;
		fillcolor=:gray,
		line_width = 0.0,
		opacity = 0.3,
	)

	add_vrect!(fig, extrema(quiet_2)...;
		fillcolor=:lightgray,
		line_width = 0.0,
		opacity = 0.3,
	)

	fig
end

# ╔═╡ 73c2366b-820d-4e51-ace4-eb04fa7e8351
signal_power = maximum(avg_power_spectrum_zoomed)

# ╔═╡ beed9cee-b8af-4e69-9407-9fa1ec398ce3
quiet_mask = let
	low1, high1 = extrema(quiet_1) .* 1e6
	low2, high2 = extrema(quiet_2) .* 1e6
	@. (low1 ≤ freqs_zoomed ≤ high1) || (low2 ≤ freqs_zoomed ≤ high2)
end;

# ╔═╡ d112613a-078c-4d4d-b856-7a5dddfe2bbe
noise_floor = mean(avg_power_spectrum_zoomed[quiet_mask])

# ╔═╡ d53d8a2d-6263-4f54-ae9f-15dcdce50f09
snr = signal_power - noise_floor

# ╔═╡ 165523d5-fe0f-437f-9dd1-1f650d1b27c3
md"""
We have an estimated SNR of **$(round(snr; digits=2)) dB**. See below for the computed values:
"""

# ╔═╡ 896b3dec-7410-4779-b04f-4457a43d61bc
md"""
# 🔧 Notebook setup
"""

# ╔═╡ 2915bcb1-cc4d-4590-8487-59693d6135ac
TableOfContents()

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Chain = "8be319e6-bccf-4806-a6f7-6fae938471bc"
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
Peaks = "18e31ff7-3703-566c-8e60-38913d67486b"
PlutoPlotly = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
Chain = "~0.6.0"
FFTW = "~1.9.0"
Peaks = "~0.5.3"
PlutoPlotly = "~0.6.3"
PlutoUI = "~0.7.68"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.5"
manifest_format = "2.0"
project_hash = "80195906885b6d13810fee17e73223f57ae16f62"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

    [deps.AbstractFFTs.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.Chain]]
git-tree-sha1 = "9ae9be75ad8ad9d26395bf625dea9beac6d519f1"
uuid = "8be319e6-bccf-4806-a6f7-6fae938471bc"
version = "0.6.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "a656525c8b46aa6a1c76891552ed5381bb32ae7b"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.30.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.ColorVectorSpace.weakdeps]
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "362a287c3aa50601b0bc359053d5c2468f0e7ce0"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.11"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "797762812ed063b9b94f6cc7742bc8883bb5e69e"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.9.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6d6219a004b8cf1e0b4dbe27a2860b8e04eba0be"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.11+0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.HashArrayMappedTries]]
git-tree-sha1 = "2eaa69a7cab70a52b9687c8bf950a5a93ec895ae"
uuid = "076d061b-32b6-4027-95e0-9a2c6f6d7e74"
version = "0.2.0"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "0f14a5456bdc6b9731a5682f439a672750a09e48"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.0.4+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "a007feb38b422fbdab534406aeca1b86823cb4d6"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "5de60bc6cb3899cd318d80d627560fae2e2d99ae"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.0.1+1"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Peaks]]
deps = ["RecipesBase", "SIMD"]
git-tree-sha1 = "75d0ce1c30696d77bc60840222d7fc5d549ebf5f"
uuid = "18e31ff7-3703-566c-8e60-38913d67486b"
version = "0.5.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PlotlyBase]]
deps = ["ColorSchemes", "Colors", "Dates", "DelimitedFiles", "DocStringExtensions", "JSON", "LaTeXStrings", "Logging", "Parameters", "Pkg", "REPL", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "28278bb0053da0fd73537be94afd1682cc5a0a83"
uuid = "a03496cd-edff-5a9b-9e67-9cda94a718b5"
version = "0.8.21"

    [deps.PlotlyBase.extensions]
    DataFramesExt = "DataFrames"
    DistributionsExt = "Distributions"
    IJuliaExt = "IJulia"
    JSON3Ext = "JSON3"

    [deps.PlotlyBase.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    JSON3 = "0f8b85d8-7281-11e9-16c2-39a750bddbf1"

[[deps.PlutoPlotly]]
deps = ["AbstractPlutoDingetjes", "Artifacts", "ColorSchemes", "Colors", "Dates", "Downloads", "HypertextLiteral", "InteractiveUtils", "LaTeXStrings", "Markdown", "Pkg", "PlotlyBase", "PrecompileTools", "Reexport", "ScopedValues", "Scratch", "TOML"]
git-tree-sha1 = "4fb7c9595eaad32d817cac8c5fa1f90daa83aa4c"
uuid = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
version = "0.6.3"

    [deps.PlutoPlotly.extensions]
    PlotlyKaleidoExt = "PlotlyKaleido"
    UnitfulExt = "Unitful"

    [deps.PlutoPlotly.weakdeps]
    PlotlyKaleido = "f2990250-8cf9-495f-b13a-cce12b45703c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "ec9e63bd098c50e4ad28e7cb95ca7a4860603298"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.68"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "fea870727142270bdf7624ad675901a1ee3b4c87"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.7.1"

[[deps.ScopedValues]]
deps = ["HashArrayMappedTries", "Logging"]
git-tree-sha1 = "1147f140b4c8ddab224c94efa9569fc23d63ab44"
uuid = "7e506255-f358-4e82-b7e4-beb19740aa63"
version = "1.3.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

    [deps.Statistics.weakdeps]
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.Tricks]]
git-tree-sha1 = "6cae795a5a9313bbb4f60683f7263318fc7d1505"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.10"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d5a767a3bb77135a99e433afe0eb14cd7f6914c3"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2022.0.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╟─ca1c0637-52c7-4168-9138-f5b91f9578dc
# ╟─525b7253-59d6-4253-a61f-75db1b61dd34
# ╠═2f45818d-34bc-4bd2-a4ac-c9552a7b144e
# ╠═158c0826-8eb4-4b0e-9d2b-65b6e02ebb44
# ╟─aeef9d9a-a111-42e6-94af-b47cc9c0c667
# ╠═b2c12094-bd4f-4fd8-9c79-2f2a351605c1
# ╠═282076e3-d070-438b-b30e-8474fa97a3b9
# ╟─762d18d8-28fa-445c-b615-87189d7e5a85
# ╟─b2aa6bf1-9943-47db-bc24-63ba3f90446e
# ╠═9b8715fb-8921-4e9f-8f46-7d68a5a69e98
# ╠═f1916076-f376-420d-8a91-00039badec1c
# ╠═6c06e8a7-9565-4ed1-8369-ccbeb240e820
# ╠═be7a8621-a11b-4b9b-b87f-ce15ea1d9236
# ╠═369e9413-d489-4a23-826c-9b14387924dc
# ╟─edb8cf40-3105-4ce6-b312-37461fb1cd23
# ╟─43e78d19-e465-43c3-9b7a-7e7e8c8c382b
# ╟─01ae457d-0093-4e92-bab4-17422ce3b97a
# ╟─0486c276-88fd-4dce-b991-e12ec75b699a
# ╠═b41854aa-29ea-4806-a12e-a3585a446539
# ╠═8bd511ab-d7af-4f68-9a29-9aaed83b3774
# ╠═cdbf8282-69b4-4a9d-9f58-e7aa270bf868
# ╟─78116562-b23d-4242-b4b9-f3bb669166f0
# ╟─78c6aebd-25b6-4890-bc6c-e8e4f63eef67
# ╟─3f4eaf24-3855-4f24-8024-b17a652aa1ad
# ╠═7e6b20a7-2a8d-48c4-8422-d371b81da6a5
# ╠═960903a9-2de8-4dc7-8883-da968562889e
# ╟─c93a0c97-f322-41ae-9f6d-ce5647b46b2d
# ╟─f5343371-806b-44d6-aa0d-f83444217ecf
# ╟─af6ddaff-6be3-44fa-bfc6-03e7ea1b4c32
# ╟─165523d5-fe0f-437f-9dd1-1f650d1b27c3
# ╠═d112613a-078c-4d4d-b856-7a5dddfe2bbe
# ╠═73c2366b-820d-4e51-ace4-eb04fa7e8351
# ╠═d53d8a2d-6263-4f54-ae9f-15dcdce50f09
# ╠═beed9cee-b8af-4e69-9407-9fa1ec398ce3
# ╟─896b3dec-7410-4779-b04f-4457a43d61bc
# ╠═2915bcb1-cc4d-4590-8487-59693d6135ac
# ╠═a7566d2f-b355-40aa-ada7-7a8830e372b3
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
