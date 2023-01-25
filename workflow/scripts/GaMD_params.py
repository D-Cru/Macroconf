# Insert GaMD parameters from sample config specification

sample_info = snakemake.wildcards.index
namespace = snakemake.params.sample
# print(snakemake.params.nstlim)

# Set total simulation time to production time given in nstlim
# + equilibration time
equil_time = 52  # time in ns
equil_steps = equil_time * 1000 / float(namespace["dt"])
prod_steps = int(namespace["nstlim"])
namespace["nstlim"] = str(int(equil_steps) + prod_steps)
print(namespace)
if namespace["igamd"] == "nan":
    namespace["igamd"] = str(3)

if int(namespace["simtime"]) > 1000:
    namespace["ntpr"] = namespace["simtime"]
    namespace["ntwx"] = namespace["simtime"]
    namespace["ntwr"] = namespace["simtime"]
else:
    namespace["ntpr"] = 1000
    namespace["ntwx"] = 1000
    namespace["ntwr"] = 1000

with open(snakemake.input.template) as f:
    params = f.read()

params_formatted = params.format(**namespace)

with open(snakemake.output.param, "w") as f:
    f.write(params_formatted)
