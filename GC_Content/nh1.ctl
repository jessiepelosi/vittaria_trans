alphabet=DNA
genetic_code=Standard

input.sequence.format=Fasta
input.sequence.file=$(file)

input.tree.file=Vittaria.partition.txt.treefile
input.tree.format=Newick
init.brlen.method=Input 


# ------------------------------------------------------------------------------
#                           Substitution model parameters
# ------------------------------------------------------------------------------

nonhomogeneous = one_per_branch
model = T92(kappa=2, theta=0.4)
# Initial values for parameters:
# kappa (K80, T92, HKY85, F84)
# kappa1, kappa2 (TN93)
# theta (T92)
# pi1, piT, piC and piG (HKY85, F84 and TN93)

nonhomogeneous_one_per_branch.shared_parameters = T92.kappa
nonhomogeneous.root_freq=GC(theta=0.4)

# Rate Across Sites variation:
# gamma or constant
rate_distribution = Invariant(dist=Gamma(n=4), p=0.2)

likelihood.recursion = simple
likelihood.recursion_simple.compression = recursive

# ----------------------------------------------------------------------------------------
#                                     Optimization
# ----------------------------------------------------------------------------------------

optimization=FullD(derivatives=Gradient)
optimization.final=powell
optimization.verbose = 1
optimization.max_number_f_eval = 100000
optimization.tolerance = 0.0000001
optimization.message_handler = ./$(file).nh1.messages
optimization.profiler = ./$(file).nh1.profile

optimization.topology = no
optimization.topology.nstep=4
optimization.topology.numfirst=no
optimization.topology.tolerance.before=100
optimization.topology.tolerance.during=100
optimization.method=fullD
optimization.method_DB.nstep=15
optimization.scale_first=no
optimization.verbose=3

# Should we write the resulting tree? none or file name.
output.tree.file = ./$(file).total.nh1.dnd
#output.tree_ids.file = ./$(file).total.nh1.ids

# Alignment information log file (site specific rates, etc):
output.infos = ./$(file).total.nh1.infos

# Write numerical parameter estimated values:
output.estimates = ./$(file).total.nh1.params.txt

output.estimates.alias=0

# ----------------------------------------------------------------------------------------
#                                     Bootstrap
# ----------------------------------------------------------------------------------------

bootstrap.number = 0
# Tell if numerical parameters should be kept to their initial value when bootstrapping: 
bootstrap.approximate = no
# Set this to yes for detailed output when bootstrapping. 
bootstrap.verbose = no
bootstrap.output.file = 
