import os
import shutil

if __name__ == '__main__':
    
    MAX_STEP = 18000000
    STEP = 18000000
    counter = 0
    FOLDER = "chi_2_0"
    HOME_DIR = os.getcwd()
    os.chdir("BD")
    os.system("make")
    os.chdir(HOME_DIR)
  
    
    steps = [100000,
             1000000, 
             2000000, 
             3000000, 
             4000000, 
             5000000, 
             6000000, 
             7000000, 
             8000000, 
             9000000, 
             10000000, 
             11000000, 
             12000000, 
             13000000, 
             14000000,
             15000000]
    distance = [21, 
                30, 
                40, 
                50, 
                60, 
                70, 
                80, 
                90, 
                100, 
                110, 
                120, 
                130, 
                140, 
                150, 
                160, 
                170]
    names = dict()
    for s, d in zip(steps, distance):
        names[s] = d
    
    with open(FOLDER+"/TRACK", 'r') as f:
        if f:
            num_atom = int(f.readline().split()[1])
        while f:
            num_step = int(f.readline().split()[1])
            # print(num_step)
            if num_step in steps:
                os.mkdir(f"{names[num_step]}")
                shutil.copy("BrownianDynamic", f"{names[num_step]}")
                shutil.copy("BONDS", f"{names[num_step]}")
                shutil.copy("CONTR", f"{names[num_step]}")
                shutil.copy("BrownianDynamic", f"{names[num_step]}")
                shutil.copy(FOLDER+"/FORCE", f"{names[num_step]}")
                coord = open(f"{names[num_step]}/COORD", 'w')
                group = open(f"{names[num_step]}/GROUP", 'w')
                coord.writelines('num_atom '+str(num_atom)+'\n')
                counter = 0
            for i in range(num_atom):
                st = f.readline()
                if num_step in steps:
                    counter += 1
                    coord.writelines(st)
                    if counter == 1:
                        group.writelines(f"GROUP 1 {' '.join(st.split()[1:])} \n")
                        group.writelines("1 \n")
                    if counter == 201:
                        group.writelines(f"GROUP 1 {' '.join(st.split()[1:])} \n")
                        group.writelines("201 \n")
            if num_step in steps:
                coord.close()
                group.close()
            if num_step >= MAX_STEP:
                break