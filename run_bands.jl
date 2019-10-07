
nprocs = 64
addprocs(nprocs - 1)

@everywhere include("/global/cscratch1/sd/d_abram/tb/bands_one_step.jl")
@everywhere using Bands


hop_file = "T300K_n6_ortho_hopping.txt"
dat = readdlm(open(hop_file),Float64)

nx = 8
ncells = nx^3
nt = Int64(round(size(dat)[1]/(ncells*10)))
@everywhere using DistributedArrays

bg = gap
gaps = @DArray [bg(it,nx,dat) for it in 1:nt]

stats = [mean(gaps), sqrt(var(gaps))]
println(stats)

writedlm("gaps.dat",stats)
writedlm("allgaps.dat", gaps)

