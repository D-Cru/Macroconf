# Insert cMD parameters from sample config specification

sample_info = snakemake.wildcards.index
namespace = snakemake.params.sample

# For simulations longer than 1000 ns, adjust the output frequency
# to output less frequently
if int(namespace["simtime"]) > 1000:
    namespace["ntpr"] = namespace["simtime"]
    namespace["ntwx"] = namespace["simtime"]
    namespace["ntwr"] = namespace["simtime"]
else:
    namespace["ntpr"] = 1000
    namespace["ntwx"] = 1000
    namespace["ntwr"] = 1000

# Read in the template
with open(snakemake.input.template) as f:
    params = f.read()

# Replace the template with the actual values
params_formatted = params.format(**namespace)

# Write the formatted template to file
with open(snakemake.output.param, "w") as f:
    f.write(params_formatted)
