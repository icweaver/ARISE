### A Pluto.jl notebook ###
# v0.20.13

#> custom_attrs = ["hide-enabled"]

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
using Chain, DataDeps, FFTW, Peaks, Statistics

# Data viz
using PlutoUI, PlutoPlotly

# Use DataDeps.jl for dataset management
# Auto-download data to current directory by default
ENV["DATADEPS_ALWAYS_ACCEPT"] = "true"
ENV["DATADEPS_LOAD_PATH"] = @__DIR__
DataDep(
	"data",
	"""
	ARISE Data Files
	Website: https://agiseti.com/
	""",
	["https://www.dropbox.com/scl/fi/6s0p463op00f0vpykpvfj/Voyager_output.ci8?rlkey=t7l2c32w91h46n9p4rj3goife&e=1&st=re8jzzpu&dl=1"],
	["ba5c565d40fef2ccc39c862bc8c0797272563c6d01eeeafc290798f86db244e6"]
) |> register
end

# ╔═╡ 0a12f6ad-26b5-4d64-ab27-844385085e9b
md"""
_This notebook is modified from the [ARISE materials](https://agiseti.com/page14.html) prepared by J. Earwicker_
"""

# ╔═╡ ca1c0637-52c7-4168-9138-f5b91f9578dc
md"""
# 🛰️ Data Science in Radio Astronomy II - Detecting Voyager 1

In this part of the lab, you will:

1. Load Voyager 1 data and compute its power spectrum.
1. Zoom in on a specified frequency range (Region of Interest, ROI).
1. Perform peak detection in the ROI.
1. Estimate the noise floor and compute the Signal-to-Noise Ratio (SNR).


!!! note "Navigating this notebook"
	* Use the table of contents located on the right-hand side to jump to different sections.
	* This notebook can be opened in multiple tabs in your browser to view and edit separate parts at the same time.
	* Clicking the "eye" icon that appears on hover at the top-left of each cell to reveal its source code. Note, this only works on local versions of the notebook. Luckily, package installation and environment setup are all nicely automated once Julia is installed on your computer. [See here](https://plutojl.org/en/docs/install/) to get started.
	* [Pluto.jl notebooks](https://plutojl.org/) are written in plain Julia; hacking is encouraged! [Source code here](https://github.com/icweaver/ARISE/blob/main/voyager/notebook.jl).
"""

# ╔═╡ 525b7253-59d6-4253-a61f-75db1b61dd34
md"""
## Step 1: Data loading 📦

1. Select either `Download dataset` (default) or `Load local dataset`
1. Verify the number of samples loaded.
1. In a separate markdown cell, note what you expect from the data (e.g., file size, data type, etc.).

!!! tip
	This dataset is fairly large (~ 300 MB), so the initial download may take a while on slower connections. This notebook will handle this automatically for you, but you can also load the file locally by selecting the `Load local dataset` option from the dropdown menu below.
"""

# ╔═╡ 542ae572-4506-41ef-afdf-e93a964e7f59
@bind choice Select(["Download dataset", "Load local dataset"];
	default = "Download dataset",
)

# ╔═╡ 86ff678e-e11c-494d-a1ae-d420a16b48e6
md"""
**Dataset location displayed below:**
"""

# ╔═╡ 7e639dd0-e9a4-43dc-be41-aeb70077ae01
begin
if choice == "Download dataset"
	# OS-agnostic file path that is automatically generated
	dataset = @datadep_str joinpath("data", "Voyager_output.ci8")
else
	# Choose pre-downladed file manually
	@bind dataset FilePicker()
end
end

# ╔═╡ d3f2343f-4ca1-4f73-bca3-450b58aff7a8
md"""
!!! note "Notes"
	👉 Your notes here
"""

# ╔═╡ deac97f2-c193-4a34-b040-ade9857fa539
load_data(dataset::AbstractString) = reinterpret(ComplexF32, read(dataset))

# ╔═╡ 335c0292-a78b-401d-b049-0139861b6374
load_data(dataset::AbstractDict) = reinterpret(ComplexF32, dataset["data"])

# ╔═╡ 2844080f-52a5-4921-8569-6f9023563b73
filtered_data = load_data(dataset)

# ╔═╡ aeef9d9a-a111-42e6-94af-b47cc9c0c667
md"""
## Step 2: Frequency analysis 💻

1. Note the range of frequencies computed in the given FFT segment.
1. In a markdown cell, write down your observations regarding the frequency range.
"""

# ╔═╡ b2c12094-bd4f-4fd8-9c79-2f2a351605c1
begin
	# Acquisition parameters
	const center_freq = 8431e6    # Center frequency in Hz
	const fs = 2.92969e6          # Sampling frequency in Hz
	const FFT_size = 16384        # Number of samples per FFT segment
end;

# ╔═╡ 282076e3-d070-438b-b30e-8474fa97a3b9
# Compute frequency axis for one FFT segment
freqs_corrected = fftshift(fftfreq(FFT_size, fs)) .+ center_freq

# ╔═╡ a7151508-50fb-4bf9-91a8-ad0497735db4
@debug "Frequency axis computed =]"

# ╔═╡ 762d18d8-28fa-445c-b615-87189d7e5a85
md"""
**Frequency range:**
$(round(first(freqs_corrected)/1e6; digits=2)) Hz to
$(round(last(freqs_corrected)/1e6; digits=2)) Hz
"""

# ╔═╡ 618ea634-949d-4851-8e0a-4154c6ab9a8c
md"""
!!! note "Notes"
	👉 Your notes here
"""

# ╔═╡ b2aa6bf1-9943-47db-bc24-63ba3f90446e
md"""
## Step 3: Entire power spectrum 📈

1. Analyze the cells used below to segment the data, compute the FFT power spectra, and plot the entire power spectrum.
1. Observe the overall spectrum.
1. In a markdown cell, describe any features you notice.
"""

# ╔═╡ e88c5f9e-cbbd-4117-a313-3687704ca728
md"""
### Segmentation

This portion segments the data into FFT-sized chunks
"""

# ╔═╡ 9b8715fb-8921-4e9f-8f46-7d68a5a69e98
num_segments = length(filtered_data) ÷ FFT_size

# ╔═╡ f1916076-f376-420d-8a91-00039badec1c
trimmed_data = filtered_data[begin : num_segments*FFT_size] # Trim extra samples

# ╔═╡ 6c06e8a7-9565-4ed1-8369-ccbeb240e820
data_segments = reshape(trimmed_data, FFT_size, num_segments)

# ╔═╡ a810ed37-ff0d-4e0c-b30c-957916768f31
md"""
### Compute FFT

We next compute FFT for each segment and obtain power spectra.
"""

# ╔═╡ be7a8621-a11b-4b9b-b87f-ce15ea1d9236
power_spectra = abs.(fftshift(fft(data_segments, 1))).^2

# ╔═╡ 5e3b3029-37e1-48a2-b653-720ba720a677
md"""
### Time average

We then compute the time-averaged power spectrum (in dB).
"""

# ╔═╡ 369e9413-d489-4a23-826c-9b14387924dc
avg_power_spectrum = 10 * log10.(mean(power_spectra, dims=2) .+ 1e-12) |> vec

# ╔═╡ 5a5a0bcd-0031-40cf-8a82-a99ead51f60a
md"""
This is a convenient way to convert the column matrix returned by

```julia
mean(power_spectra, dims=2)
```

into a simple, plotly-friendly vector.
""" |> x -> details(md"What is `vec` used for?", x)

# ╔═╡ 6a9d9687-e897-4e92-8bba-b577f2c3b5da
md"""
### Plot

Finally, we plot the entire power spectrum below.
"""

# ╔═╡ edb8cf40-3105-4ce6-b312-37461fb1cd23
let
	p = scatter(x=freqs_corrected/1e6, y=avg_power_spectrum)

	layout = Layout(
		xaxis = attr(title="Frequency (MHz)"),
		yaxis = attr(title="Power (dB)"),
		title = "Voyager 1 Entire Power Spectrum",
	)

	
	plot(p, layout)
end

# ╔═╡ aa5666a4-432a-4f94-aed9-fd5c25ba2a1c
md"""
!!! note "Notes"
	👉 Your notes here
"""

# ╔═╡ 43e78d19-e465-43c3-9b7a-7e7e8c8c382b
md"""
## Step 4: Region of Interest (ROI) 🔎

1. Modify the ROI boundaries in the input boxes below the plot if desired. The relevant computation cells will update automatically.
1. Observe the zoomed-in spectrum.
1. Write a short note (in a markdown cell) about the features visible in the ROI.
"""

# ╔═╡ 0486c276-88fd-4dce-b991-e12ec75b699a
md"""
**Lower bound:** $(@bind roi_min NumberField(8429:0.1:8431; default=8430.1)) MHz,
**Upper bound:** $(@bind roi_max NumberField(8431:0.1:8433; default=8431.9)) MHz
"""

# ╔═╡ dce75344-52f8-4681-96c9-56bf25fe2d9e
md"""
Below are the relevant computations made to produce the above plot.
"""

# ╔═╡ b41854aa-29ea-4806-a12e-a3585a446539
# Create a mask for the ROI using the frequency axis computed earlier
freq_mask = @. roi_min*1e6 ≤ freqs_corrected ≤ roi_max*1e6

# ╔═╡ 87df2143-b66f-453b-a8db-6b386ba9f796
md"""
This is a convenient way to broadcast functions in Julia. For example, the following two lines provided equivalent results:

```julia
[6, 7, 8].^2 .+ 1 # [37, 50, 65]
@. [6, 7, 8]^2 + 1 # [37, 50, 65]
```
	
[See here](https://docs.julialang.org/en/v1/manual/functions/#man-vectorized) for more.
""" |> x -> details(md"What's `@.`?", x)

# ╔═╡ 8bd511ab-d7af-4f68-9a29-9aaed83b3774
freqs_zoomed = @views freqs_corrected[freq_mask]

# ╔═╡ 829aa3d7-3d84-4837-9b5b-846532ad0825
md"""
These are nice ways to "view" into portions of the data, instead of creating a separate copy in memory each time, which can be more computationally expensive. [See here](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-views) for more.
""" |> x -> details(md"What are `views`?", x)

# ╔═╡ cdbf8282-69b4-4a9d-9f58-e7aa270bf868
# Compute the time-averaged power spectrum for the ROI (in dB)
avg_power_spectrum_zoomed = @views 10 * log10.(mean(power_spectra[freq_mask, :], dims=2) .+ 1e-12) |> vec

# ╔═╡ 01ae457d-0093-4e92-bab4-17422ce3b97a
let
	p = scatter(x=freqs_zoomed/1e6, y=avg_power_spectrum_zoomed)

	layout = Layout(
		xaxis = attr(title="Frequency (MHz)"),
		yaxis = attr(title="Power (dB)"),
		title = "Voyager 1 Power Spectrum (ROI)",
	)

	plot(p, layout)
end

# ╔═╡ cbcaca71-c2db-4227-87f0-a387b837d5f6
md"""
!!! note "Notes"
	👉 Your notes here
"""

# ╔═╡ 78116562-b23d-4242-b4b9-f3bb669166f0
md"""
## Step 5: Peak detection 🗻

1. Adjust the peak detection parameters in the table below if necessary.
1. Verify the detected peaks visually on the zoomed in plot.
1. In a markdown cell, explain why certain peaks were detected and how the parameters affect the results.
"""

# ╔═╡ 3f4eaf24-3855-4f24-8024-b17a652aa1ad
md"""
param | value | description
:--   | :--   | :--
peak\_height\_offset | $(@bind peak_height_offset NumberField(1:0.1:100; default=30)) | dB offset below the ROI's maximum power to set threshold
min_distance | $(@bind min_distance NumberField(1:100; default=10)) | Minimum separation between peaks (in index units)
prominence | $(@bind prominence NumberField(1:0.1:100; default=5)) | Minimum prominence (in dB) for a peak
"""

# ╔═╡ 10f88e0b-7764-486b-ba95-40b158b4c2fd
md"""
Below are the relevant computations made to produce the above plot.
"""

# ╔═╡ 960903a9-2de8-4dc7-8883-da968562889e
# Calculate the detection threshold based on the maximum value in the ROI spectrum
peak_height = maximum(avg_power_spectrum_zoomed) - peak_height_offset

# ╔═╡ 7e6b20a7-2a8d-48c4-8422-d371b81da6a5
# Detect peaks within the ROI spectrum
pks = @chain avg_power_spectrum_zoomed begin
	findmaxima(min_distance)
	peakproms(min=prominence)
	peakheights(min=peak_height)
end

# ╔═╡ 78c6aebd-25b6-4890-bc6c-e8e4f63eef67
let
	p1 = scatter(x=freqs_zoomed/1e6, y=avg_power_spectrum_zoomed;
		name = "Spectrum ROI",
	)

	mask = pks.indices
	p2 = scatter(x=freqs_zoomed[mask]/1e6, y=avg_power_spectrum_zoomed[mask];
		mode = :markers,
		marker = attr(symbol=:x),
		name = "Detected Peaks",
	)

	layout = Layout(
		xaxis = attr(title="Frequency (MHz)"),
		yaxis = attr(title="Power (dB)"),
		title = "Peak Detection on ROI",
	)
	
	plot([p1, p2], layout)
end

# ╔═╡ 5e1590ab-e2b3-48aa-b2a8-9514b20e3dc7
md"""
A nice way to chain functions together, especially when they start to get long. The following would produce equivalent results:

```julia
floor(exp(log10(sqrt(2)))) # 1.0
	
2 |> sqrt |> log10 |> exp |> floor # 1.0

@chain 2 begin
	sqrt
	log10
	exp
	floor
end # 1.0
```

[See here](https://github.com/jkrumbiegel/Chain.jl) for more.
""" |> x -> details(md"What's `@chain`?", x)

# ╔═╡ 3bb02848-4580-4644-a2b8-87c05cf3e60e
md"""
!!! note "Notes"
	👉 Your notes here
"""

# ╔═╡ c93a0c97-f322-41ae-9f6d-ce5647b46b2d
md"""
## Step 6: SNR computation 📢

1. Using the sliders below, define quiet regions where no signal is expected.
1. Estimate the noise floor by averaging the power in these quiet regions.
1. Compute the Signal-to-Noise Ratio (SNR) as the difference between the main signal power and the noise floor.
1. In a markdown cell, discuss alternative methods (e.g., using the median) for noise estimation.
"""

# ╔═╡ f5343371-806b-44d6-aa0d-f83444217ecf
md"""
**Quiet region 1:** $(@bind quiet_1 RangeSlider(8430:0.1:8431; default=8430.4:0.1:8430.9)) MHz

**Quiet region 2:** $(@bind quiet_2 RangeSlider(8431:0.1:8432; default=8431.2:0.1:8431.6)) MHz
"""

# ╔═╡ af6ddaff-6be3-44fa-bfc6-03e7ea1b4c32
let
	p1 = scatter(x=freqs_zoomed/1e6, y=avg_power_spectrum_zoomed;
		name = "Spectrum ROI",
	)

	
	layout = Layout(
		xaxis = attr(title="Frequency (MHz)"),
		yaxis = attr(title="Power (dB)"),
		title = "Quiet Regions for Noise Floor Estimation",
	)

	fig = plot(p1, layout)

	add_vrect!(fig, extrema(quiet_1)...;
		fillcolor=:gray,
		line_width = 0.0,
		opacity = 0.3,
		showlegend = true,
		name = "Quiet Region 1"
	)

	add_vrect!(fig, extrema(quiet_2)...;
		fillcolor=:lightgray,
		line_width = 0.0,
		opacity = 0.3,
		showlegend = true,
		name = "Quiet Region 2",
	)

	fig
end

# ╔═╡ beed9cee-b8af-4e69-9407-9fa1ec398ce3
# Create a mask for the quiet regions
quiet_mask = let
	low1, high1 = extrema(quiet_1) .* 1e6
	low2, high2 = extrema(quiet_2) .* 1e6
	@. (low1 ≤ freqs_zoomed ≤ high1) || (low2 ≤ freqs_zoomed ≤ high2)
end;

# ╔═╡ d112613a-078c-4d4d-b856-7a5dddfe2bbe
# Estimate the noise floor by averaging the power in the quiet regions (in dB)
noise_floor = mean(avg_power_spectrum_zoomed[quiet_mask])

# ╔═╡ 73c2366b-820d-4e51-ace4-eb04fa7e8351
# Identify the main signal peak in the ROI (assumed to be the maximum power)
signal_power = maximum(avg_power_spectrum_zoomed)

# ╔═╡ d53d8a2d-6263-4f54-ae9f-15dcdce50f09
# Compute SNR (in dB)
snr = signal_power - noise_floor

# ╔═╡ 165523d5-fe0f-437f-9dd1-1f650d1b27c3
md"""
We have an estimated SNR of **$(round(snr; digits=2)) dB**. See below for the relevant computations:
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
DataDeps = "124859b0-ceae-595e-8997-d05f6a7a8dfe"
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
Peaks = "18e31ff7-3703-566c-8e60-38913d67486b"
PlutoPlotly = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
Chain = "~0.6.0"
DataDeps = "~0.7.13"
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
project_hash = "196cd2e3038b77c63a5ca25982b8ae8f359792ae"

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

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.Chain]]
git-tree-sha1 = "9ae9be75ad8ad9d26395bf625dea9beac6d519f1"
uuid = "8be319e6-bccf-4806-a6f7-6fae938471bc"
version = "0.6.0"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

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

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "d9d26935a0bcffc87d2613ce14c527c99fc543fd"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.5.0"

[[deps.DataDeps]]
deps = ["HTTP", "Libdl", "Reexport", "SHA", "Scratch", "p7zip_jll"]
git-tree-sha1 = "8ae085b71c462c2cb1cfedcb10c3c877ec6cf03f"
uuid = "124859b0-ceae-595e-8997-d05f6a7a8dfe"
version = "0.7.13"

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

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "d36f682e590a83d63d1c7dbd287573764682d12a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.11"

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

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "PrecompileTools", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "ed5e9c58612c4e081aecdb6e1a479e18462e041e"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.17"

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

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "f02b56007b064fbfddb4c9cd60161b6dd0f40df3"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.1.0"

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

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

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

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "f1a7e086c677df53e064e0fdd2c9d0b0833e3f6e"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.5.0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "87510f7292a2b21aeff97912b0898f9553cc5c2c"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.1+0"

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

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

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

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

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
# ╟─0a12f6ad-26b5-4d64-ab27-844385085e9b
# ╟─ca1c0637-52c7-4168-9138-f5b91f9578dc
# ╟─525b7253-59d6-4253-a61f-75db1b61dd34
# ╟─542ae572-4506-41ef-afdf-e93a964e7f59
# ╟─86ff678e-e11c-494d-a1ae-d420a16b48e6
# ╟─7e639dd0-e9a4-43dc-be41-aeb70077ae01
# ╟─2844080f-52a5-4921-8569-6f9023563b73
# ╠═d3f2343f-4ca1-4f73-bca3-450b58aff7a8
# ╟─deac97f2-c193-4a34-b040-ade9857fa539
# ╟─335c0292-a78b-401d-b049-0139861b6374
# ╟─aeef9d9a-a111-42e6-94af-b47cc9c0c667
# ╠═b2c12094-bd4f-4fd8-9c79-2f2a351605c1
# ╠═282076e3-d070-438b-b30e-8474fa97a3b9
# ╟─a7151508-50fb-4bf9-91a8-ad0497735db4
# ╟─762d18d8-28fa-445c-b615-87189d7e5a85
# ╠═618ea634-949d-4851-8e0a-4154c6ab9a8c
# ╟─b2aa6bf1-9943-47db-bc24-63ba3f90446e
# ╟─e88c5f9e-cbbd-4117-a313-3687704ca728
# ╠═9b8715fb-8921-4e9f-8f46-7d68a5a69e98
# ╠═f1916076-f376-420d-8a91-00039badec1c
# ╠═6c06e8a7-9565-4ed1-8369-ccbeb240e820
# ╟─a810ed37-ff0d-4e0c-b30c-957916768f31
# ╠═be7a8621-a11b-4b9b-b87f-ce15ea1d9236
# ╟─5e3b3029-37e1-48a2-b653-720ba720a677
# ╠═369e9413-d489-4a23-826c-9b14387924dc
# ╟─5a5a0bcd-0031-40cf-8a82-a99ead51f60a
# ╟─6a9d9687-e897-4e92-8bba-b577f2c3b5da
# ╟─edb8cf40-3105-4ce6-b312-37461fb1cd23
# ╠═aa5666a4-432a-4f94-aed9-fd5c25ba2a1c
# ╟─43e78d19-e465-43c3-9b7a-7e7e8c8c382b
# ╟─01ae457d-0093-4e92-bab4-17422ce3b97a
# ╟─0486c276-88fd-4dce-b991-e12ec75b699a
# ╟─dce75344-52f8-4681-96c9-56bf25fe2d9e
# ╠═b41854aa-29ea-4806-a12e-a3585a446539
# ╟─87df2143-b66f-453b-a8db-6b386ba9f796
# ╠═8bd511ab-d7af-4f68-9a29-9aaed83b3774
# ╟─829aa3d7-3d84-4837-9b5b-846532ad0825
# ╠═cdbf8282-69b4-4a9d-9f58-e7aa270bf868
# ╠═cbcaca71-c2db-4227-87f0-a387b837d5f6
# ╟─78116562-b23d-4242-b4b9-f3bb669166f0
# ╟─3f4eaf24-3855-4f24-8024-b17a652aa1ad
# ╟─78c6aebd-25b6-4890-bc6c-e8e4f63eef67
# ╟─10f88e0b-7764-486b-ba95-40b158b4c2fd
# ╠═960903a9-2de8-4dc7-8883-da968562889e
# ╠═7e6b20a7-2a8d-48c4-8422-d371b81da6a5
# ╟─5e1590ab-e2b3-48aa-b2a8-9514b20e3dc7
# ╠═3bb02848-4580-4644-a2b8-87c05cf3e60e
# ╟─c93a0c97-f322-41ae-9f6d-ce5647b46b2d
# ╟─f5343371-806b-44d6-aa0d-f83444217ecf
# ╟─af6ddaff-6be3-44fa-bfc6-03e7ea1b4c32
# ╟─165523d5-fe0f-437f-9dd1-1f650d1b27c3
# ╠═beed9cee-b8af-4e69-9407-9fa1ec398ce3
# ╠═d112613a-078c-4d4d-b856-7a5dddfe2bbe
# ╠═73c2366b-820d-4e51-ace4-eb04fa7e8351
# ╠═d53d8a2d-6263-4f54-ae9f-15dcdce50f09
# ╟─896b3dec-7410-4779-b04f-4457a43d61bc
# ╠═2915bcb1-cc4d-4590-8487-59693d6135ac
# ╠═a7566d2f-b355-40aa-ada7-7a8830e372b3
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
