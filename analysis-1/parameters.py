import os
import sys

def getParameters(osargs_list):
	parameters = {}

	# define specific os arguments if provided
	for osarg in osargs_list:
		try: parameters[osarg]=os.environ[osarg]
		except: print(f'{osarg} not defined in os')

	# iterate over sys.args and override any provided arguments
	for sysarg in sys.argv[1:]:
		name, value=sysarg.split('=')
		parameters[name]=value
	
	return(parameters)
