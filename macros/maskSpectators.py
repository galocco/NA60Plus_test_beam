import ROOT
import argparse
import re
import shutil


# it would be nice to 
#
#
#

def read_geometry(file_path):

    #example of block of stuff for an ALPIDE
    """
        [ALPIDE_0]
        coordinates = "cartesian"
        material_budget = 0.001
        number_of_pixels = 1024, 512
        orientation = 0deg,180deg,0deg
        orientation_mode = "xyz"
        pixel_pitch = 29.24um,26.88um
        position = 0um,0um,0um
        spatial_resolution = 5um,5um
        time_resolution = 2us
        type = "alpide"
    """

    z_position = []
    alpide_list = []

    try:
        with open(file_path, 'r') as file:
            for line in file:
                if "[ALPIDE_" in line:
                    alpide_list.append(float(re.search(r'\d+', input_string)))
                if "position" in line:
                    z_position.append(float(re.findall(r'(\d+\.\d+)um', input_string)[-1]))

    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

    return z_position[2:], alpide_list[2:]

def create_masking(args):
    MERGE = args.merge
    MASK_CENTER = args.center
    MASK_EDGES = args.edges
    target_thickness = args.thickness
    radius0 = args.radius
    beam_sigma_max = args.beam_sigma
    nsigma = args.nsigma

    if args.setup == 1:
        #path of the file with the hitmap
        path = "/home/giacomo/its-corryvreckan-tools/output/2023-10_SPS/maskcreation_422231250231017232859.root"
        mask_dir = "/masking/setup1_2023"
        geometry_file = "/home/giacomo/its-corryvreckan-tools/output/2023-10_SPS/2023-10_SPS_8REF_setup1.conf"
    elif args.setup == 2:
        #path of the file with the hitmap
        path = "/home/giacomo/its-corryvreckan-tools/output/2023-10_SPS/maskcreation_422231250231017232859.root"
        mask_dir = "/masking/setup2_2023"
        geometry_file = "/home/giacomo/its-corryvreckan-tools/output/2023-10_SPS/2023-10_SPS_7REF_setup2.conf"
    elif args.setup == 3:
        #path of the file with the hitmap
        path = "/home/giacomo/its-corryvreckan-tools/output/2023-10_SPS/maskcreation_422231250231017232859.root"
        mask_dir = "/masking/setup3_2023"
        geometry_file = "/home/giacomo/its-corryvreckan-tools/output/2023-10_SPS/2023-10_SPS_8REF_setup3.conf"

    z_position, alpide_list = read_geometry(geometry_file)

    z_target = z_position[0]-25

    pixel_pitch = [29.24,26.88]

    file = ROOT.TFile(path)


    target_distance_last = z_position[-1] -z_target -target_thickness/2.

    npixels_x = 1024
    npixels_y = 512
    test_file = ROOT.TFile("test.root","recreate")
    for alpide,z_alpide in zip(alpide_list,z_position):
        if MASK_CENTER:
            th2 = ROOT.TH2D("th2_"+str(alpide),";x;y",1024,-0.5,1023.5,512,-0.5,511.5)
            mode = ""
            if MERGE:
                mode = "a"
                
                source_file = '/home/giacomo/its-corryvreckan-tools/output/ALPIDE_'+str(alpide)+'mask_ALPIDE_'+str(alpide)+'.txt'

                shutil.copy(source_file, mask_dir+"/masking_"+str(alpide)+".txt")
            f = open(mask_dir+"/masking_"+str(alpide)+".txt", mode)
            print("masking_"+str(alpide)+".txt")
            hitmap = file.Get("ClusteringSpatial/ALPIDE_"+str(alpide)+"/clusterPositionLocal")

            mean_x = hitmap.GetMean(1)
            mean_y = hitmap.GetMean(2)
            print("ALPIDE:")
            print("local:")
            print("mean x: ",mean_x)
            print("mean y: ",mean_y)
            hitmap = file.Get("ClusteringSpatial/ALPIDE_"+str(alpide)+"/clusterPositionGlobal")

            mean_x = hitmap.GetMean(1)
            mean_y = hitmap.GetMean(2)
            print("global:")
            print("mean x: ",mean_x)
            print("mean y: ",mean_y)
            theta_list = []

            radius = ROOT.TMath.Sqrt(radius0*(z_alpide-z_target)/125)**2+4*beam_sigma_max**2)
            X = radius/pixel_pitch[0]
            Y = radius/pixel_pitch[1]

            for ix in range(0,npixels_x):
                for iy in range(0,npixels_y):
                    r2 = (ix-mean_x)**2/X**2+(iy-mean_y)**2/Y**2
                    if r2 < 1:
                        f.write("p	"+str(ix)+"  "+str(iy)+"\n")
                        th2.SetBinContent(ix,iy,1)
            f.close()
        #work in progress    
        if False:
            if alpide == 5:
                d1 = npixels_y - mean_y
                if d1 > mean_y:
                    d1 = mean_y
                
                X = d1*pixel_pitch[1]/pixel_pitch[0]
                Y = d1

                for ix in range(0,npixels_x):
                    for iy in range(0,npixels_y):
                        r2 = (ix-mean_x)**2/X**2+(iy-mean_y)**2/Y**2
                        if r2 > 1:
                            f.write("p	"+str(ix)+"  "+str(iy)+"\n")
                            th2.SetBinContent(ix,iy,1)
                th2.Write()

        if MASK_EDGES and alpide != alpide_list[-1]:
            if False:
                target_distance = z_alpide -z_target -target_thickness/2.
                Lxi = int(target_distance/target_distance_last * npixels_x / 2.)
                Lyi = int(target_distance/target_distance_last * npixels_y / 2.)
                print("ALPIDE:",alpide)
                print(Lxi)
                print(Lyi)
                for iy in range(0,npixels_y):
                    for ix in range(0, Lxi+1):
                        f.write("p	"+str(ix)+"  "+str(iy)+"\n")
                    for ix in range(int(npixels_x/2.)+Lxi,npixels_x):
                        f.write("p	"+str(ix)+"  "+str(iy)+"\n")

                for ix in range(Lxi+1, int(npixels_x/2.)+Lxi+1):
                    for iy in range(0,Lyi+1):
                        f.write("p	"+str(ix)+"  "+str(iy)+"\n")
                    for iy in range(int(npixels_y/2.)+Lyi, npixels_y):
                        f.write("p	"+str(ix)+"  "+str(iy)+"\n")
    test_file.Close()


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--merge', help='Merge with other masking', action='store_true')
    parser.add_argument('-c', '--center', help='Mask the center of the ALPIDEs', action='store_true')
    parser.add_argument('-e', '--edges', help='Mask the edge of the ALPIDEs', action='store_true') #not implemented yet
    parser.add_argument('-t','--thickness', type=float, help='Set target thickness in mm', default = 5)
    parser.add_argument('-r','--radius', type=float, help='Set masking r0 in mm', default = 1)
    parser.add_argument('-b','--beam_sigma', type=float, help='Set beam size in mm', default = 0.23)
    parser.add_argument('-n','--nsigma', type=float, help='Set nsigma in mm', default = 2)
    parser.add_argument('-s','--setup', type=float, choices=[1,2,3], help='Choose between option1 and option2', default=1)
    args = parser.parse_args()

    create_masking(args)
