from corpse import simpleFastcore
from cobra.io import read_sbml_model
from cobra.io import write_sbml_model
from pandas import read_csv
from os import listdir
from os.path import isfile, join
import re
from os.path import exists
from pathlib import Path
home = str(Path.home())

# path with csv of core reactions from proteomic data (see: core-rxn.r)
path = home+"csm/"
files = [f for f in listdir(path) if f.endswith(".csv")]

for f in files:
	if exists(path+"/"+re.sub("csv$","xml",f)):
		print("Already done: "+f)
		continue
	rxn_list = read_csv(path+f)["rxn.active"].to_list()	
	org = re.findall(".*?(?=-)", f)[0]
	print("Preparing model for: "+org+" "+f)
	mod = read_sbml_model(path+"../out/"+org+".xml.gz")
	mod.solver = "cplex"
	fast_mod = simpleFastcore(model = mod, core_set = rxn_list)
	core_mod = fast_mod.run()
	write_sbml_model(core_mod, path + "/"+re.sub("csv$","xml",f))
