### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ d24d4a18-0566-4d4b-b64b-464be7cb8237
struct Item
	cost::Int
	criticality::Int
end


# ╔═╡ 92ab830e-0aa2-4ef6-a4bb-25aee2ba7799
begin
	itm1 = Item(1000,5)
	itm2 = Item();
end


# ╔═╡ 826a8142-12a9-49a3-9b2b-b0d7448274c9
num_cov_levels = 3;

# ╔═╡ 3ceb7f71-079e-4672-895c-aa0f04bd80a1
x = collect(LinRange(1,num_cov_levels,num_cov_levels))

# ╔═╡ 943665a0-4a86-43a5-b193-b2d0a1682050
y = collect(LinRange(1,0,num_cov_levels))

# ╔═╡ 2589c79d-adb0-4da6-8e4f-23b96c0f4f42
cov_lvls = (x .=> y)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.3"
manifest_format = "2.0"
project_hash = "da39a3ee5e6b4b0d3255bfef95601890afd80709"

[deps]
"""

# ╔═╡ Cell order:
# ╠═d24d4a18-0566-4d4b-b64b-464be7cb8237
# ╠═92ab830e-0aa2-4ef6-a4bb-25aee2ba7799
# ╠═826a8142-12a9-49a3-9b2b-b0d7448274c9
# ╠═3ceb7f71-079e-4672-895c-aa0f04bd80a1
# ╠═943665a0-4a86-43a5-b193-b2d0a1682050
# ╠═2589c79d-adb0-4da6-8e4f-23b96c0f4f42
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
