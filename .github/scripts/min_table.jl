# Build a PR-comment markdown table from benchpkg result JSONs, reporting the
# *minimum* time per benchmark instead of AirspeedVelocity's median. The minimum
# is robust to noisy-neighbour contention on shared CI runners (it captures the
# least-disturbed sample), so it doesn't show phantom regressions when the head
# build happens to land on a busier runner window.
#
# Usage:  julia min_table.jl <input_dir> <pkg> <base_rev> <head_rev> [label]
# The optional <label> (e.g. the runner OS) is appended to the table heading so a
# multi-platform matrix can post one distinct comment per platform.
using JSON3

input_dir, pkg, base_rev, head_rev = ARGS[1], ARGS[2], ARGS[3], ARGS[4]
label = length(ARGS) >= 5 ? ARGS[5] : ""

readrev(rev) = open(joinpath(input_dir, "results_$(pkg)@$(rev).json"), "r") do io
    JSON3.read(io, Dict{String,Any})
end

# Recursively collect leaf trials (nodes carrying a "times" array) keyed by "a/b/c".
function leaves!(acc, node, prefix = "")
    if haskey(node, "times")
        acc[prefix] = node
    elseif haskey(node, "data")
        for (k, v) in node["data"]
            name = isempty(prefix) ? String(k) : prefix * "/" * String(k)
            leaves!(acc, v, name)
        end
    end
    acc
end

mintime(node) = minimum(Float64.(node["times"]))   # ns

function fmt_time(ns)
    unit, div = ns < 1e3 ? ("ns", 1.0) :
                ns < 1e6 ? ("μs", 1e3) :
                ns < 1e9 ? ("ms", 1e6) : ("s", 1e9)
    string(round(ns / div; sigdigits = 3), " ", unit)
end

fmt_mem(node) = string(round(Int, node["allocs"]), " allocs: ", round(Int, node["memory"]), " B")

base = leaves!(Dict{String,Any}(), readrev(base_rev))
head = leaves!(Dict{String,Any}(), readrev(head_rev))

# Stable ordering: sort names, push time_to_load to the end.
names = sort(collect(keys(base)))
filter!(!=("time_to_load"), names)
haskey(base, "time_to_load") && push!(names, "time_to_load")

headlbl = length(head_rev) >= 8 ? head_rev[1:8] * "…" : head_rev

io = IOBuffer()
title = isempty(label) ? "Benchmark Results (minimum time)" :
                         "Benchmark Results (minimum time) — $label"
println(io, "## $title")
println(io)
println(io, "Reporting the **minimum** over all samples (robust to shared-runner contention), ",
            "not the median. Ratio = $base_rev / $headlbl: **>1 means the PR is faster**. ",
            "✅ ≥ 5 % faster, ⚠️ ≥ 5 % slower.")
println(io)

# --- time table ---
println(io, "<details open><summary>Time benchmarks</summary>")
println(io)
println(io, "|  | $base_rev | $headlbl | $base_rev / $headlbl |")
println(io, "|:--|--:|--:|--:|")
for n in names
    haskey(head, n) || continue
    mb, mh = mintime(base[n]), mintime(head[n])
    ratio = round(mb / mh; sigdigits = 3)
    flag = ratio >= 1.05 ? " ✅" : ratio <= 0.95 ? " ⚠️" : ""
    println(io, "| $n | $(fmt_time(mb)) | $(fmt_time(mh)) | $(ratio)$(flag) |")
end
println(io)
println(io, "</details>")
println(io)

# --- memory table ---
println(io, "<details><summary>Memory benchmarks</summary>")
println(io)
println(io, "|  | $base_rev | $headlbl |")
println(io, "|:--|--:|--:|")
for n in names
    haskey(head, n) || continue
    println(io, "| $n | $(fmt_mem(base[n])) | $(fmt_mem(head[n])) |")
end
println(io)
println(io, "</details>")

print(String(take!(io)))
