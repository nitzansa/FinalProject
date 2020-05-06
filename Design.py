
from tkinter import *
#gui- create the view
from tkinter import filedialog
#from Artifacts import Artifact

import sys
import os

if os.environ.get('DISPLAY','') == '':
    print('no display found. Using :0.0')
    os.environ.__setitem__('DISPLAY', ':0.0')

root = Tk()
root.geometry('800x500')
root.title("Analysis Of Bacterial Gene Cluster")


title_lbl = Label(root, text="Analysis Of Bacterial Gene Cluster", font=("Arial Bold", 20)).place(x=150, y=50)
#title_lbl.grid(column=200, row=2)

insertP_lbl = Label(root, text="Enter bacteria genomes directory path:",font=("Arial", 10)).place(x= 20, y=150)

txt_path = StringVar(root, value='Enter directory')
path_txt = Entry(root,width=30, textvariable=txt_path).place(x=255,y=150)



def browse_button():
    filename = filedialog.askdirectory()
    txt_path.set(filename)
    print(filename)
    #lbl.configure(text="Button was clicked !!")
    ##### open dialoge to enter path
path_btn = Button(root, text="browse path", command=browse_button, font=("Arial", 10)).place(x=450,y=150)


def showClusters():
    print("CH-HIT Clusters HERE")


showClusters_btn = Button(root, text="Show CD-HIT Algorithm Output", command=showClusters, font=("Arial", 10)).place(x=200,y=220)


def showProblematicStrains():
    print("Problematic Strains HERE")
   # a = Artifact("/home/local/BGU-USERS/sabagnit/CD_HIT_output_sqeuence")
    #a.variableLength()
    #a.getGenesPerCluster()
    #a.getStrainsPerCluster()
    # a.getMinStrainsPerCluster(2)
    # print(a.getSingleClusters())
    # print(len(a.getSingleClusters()))
    # a.getStrainsPerCluster()
  #  a.calcAverageMemberPerCluster()
    # a.clustersPerCountOfStrains()
   # a.downloadReport()


showProblematicStrains_btn = Button(root, text="Show Problematic Strains", command=showProblematicStrains, font=("Arial", 10)).place(x=215,y=260)










root.mainloop()