#!/home/user/venv/bin/python3

#A utility that prints out the number of hydrogen bonds between different strands in the system 

# import base
try:
    import numpy as np
except:
    import mynumpy as np
import os.path
import sys
import readers 
import subprocess
import tempfile
import os


def make_bonds_file():

	command_for_data =  'analysis_data_output_1 = { \n name = stdout \n print_every = 1 \n col_1 = { \n type=pair_energy \n} \n}'
	# PROCESSPROGRAM = os.path.join(os.path.dirname(__file__), "../build/bin/DNAnalysis")
	PROCESSPROGRAM = os.path.join("/home/user/SA-EDS/oxDNA/build/bin", "DNAnalysis")
	# PROCESSPROGRAM = os.path.join("oxDNA/build/bin", "DNAnalysis")


	if (len(sys.argv) < 3):
		print('Usage %s input_file trajectory_file topology_file [confid]' % sys.argv[0])
		sys.exit()

	confid = 0
	#now get topology file name:
	inputfile = sys.argv[1]
	conffile = sys.argv[2]
	topologyfile = sys.argv[3]
	# print("debug", inputfile, conffile, topologyfile)
	if len(sys.argv) >= 5:
		confid = int(sys.argv[4])

	# topologyfile = ""
	fin = open(inputfile)
	# for line in fin:
	#     line = line.lstrip()
	#     if not line.startswith('#'):
	# 	    if "topology" in line:
	# 		    topologyfile = line.split('=')[1].replace(' ','').replace('\n','')
	myreader = readers.LorenzoReader(conffile,topologyfile)
	mysystem = myreader.get_system()

	if not os.path.isfile(PROCESSPROGRAM):
		print("Cannot execute output_bonds program. Please make sure to compile DNAnalysis in ../bin/ directory")
		sys.exit(1)

	counter = 0



	tempfile_obj = tempfile.NamedTemporaryFile()
	launchcommand = inputfile + ' trajectory_file='+tempfile_obj.name+' '+command_for_data

	launchargs = [PROCESSPROGRAM, inputfile ,'trajectory_file='+tempfile_obj.name,command_for_data]
	launchargs = ["/bin/bash", PROCESSPROGRAM, inputfile ,'trajectory_file='+tempfile_obj.name,command_for_data]
	launchargs = ["." + PROCESSPROGRAM, inputfile ,'trajectory_file='+tempfile_obj.name,command_for_data]
	# launchargs = [PROCESSPROGRAM,"-v", inputfile ,'trajectory_file='+tempfile_obj.name,command_for_data]
	#print command_for_data
	#launchargs = [PROCESSPROGRAM,inputfile ,'trajectory_file='+conffile,command_for_data]

	bonds_file_name = "/".join(inputfile.split("/")[:-1]) + "/" + "bonds"
	bonds_file = open(bonds_file_name, 'w')

	while mysystem != False:
		mysystem.map_nucleotides_to_strands()
		mysystem.print_lorenzo_output(tempfile_obj.name,'/dev/null')
		tempfile_obj.flush()

		# print("mysystem : ", mysystem._conf)
		# launchcommand = 'trajectory_file='+tempfile_obj.name+' '+command_for_data
		# os.system(PROCESSPROGRAM+' '+launchcommand)
		# launchargs = [PROCESSPROGRAM,launchcommand]
		if counter == confid:
			# /usr/bin/bash
			myinput = subprocess.Popen(launchargs,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			subprocess.run(["chmod", "+x", "/home/user/SA-EDS/oxDNA/build/bin/DNAnalysis"])
			mystdout,mystderr = myinput.communicate()
			str_mystdout = mystdout.decode("utf-8")
			str_mystderr = mystderr.decode("utf-8")
			# /home/user/SA-EDS/oxDNA/build/bin/DNAnalysis: /home/user/SA-EDS/oxDNA/build/bin/DNAnalysis: cannot execute binary file
			
			bonds_file.write(str_mystdout)
			sys.exit(1)
		counter += 1
		mysystem = myreader.get_system()

if __name__ == '__main__':
	make_bonds_file()
	pass
