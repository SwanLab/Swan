from utils import exec1
import os   
dirs = os.listdir('.')

for dir in [x for x in dirs if os.path.isdir(x)]:  
    files = os.listdir(dir)
    for file in [f for f in files if f.endswith('.py')]:    
        exec1("python "+dir+"/"+file)

