from ROOT import TH1D,TH2D, TFile, TF2, TCanvas, kRed, TLegend, TFile, gStyle, kRainBow, TTree
import ROOT
import uproot
import numpy as np
import math
import argparse
import vertexer
import pandas as pd
import yaml
import os
import time
from loading_bar import *
import warnings

# Ignore all RuntimeWarnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

def read_tree(path, compute_sigmas = False, config = "", nsigma = 1, suffix ="", update = False, run_vertexing = True, run_analysis = True):
    #put here the path of the file you want to analyze

    df_beam = uproot.open(path)['Tracking4D/tree/event/beam'].arrays(library="pd")
    df_tracks = uproot.open(path)['Tracking4D/tree/event/tracks'].arrays(library="pd", entry_stop=20000)
 
    with open(os.path.expandvars(config), 'r') as stream:
        try:
            params = yaml.full_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    
    target = params["TARGET"]
    z_target = params["Z_TARGET"]
    thickness = params["THICKNESS"]
    interaction_prob = params["INTERACTION_PROB"]

    sigma_x = params["SIGMA_X_Pb"]
    sigma_y = params["SIGMA_Y_Pb"]
    mean_x = params["MEAN_X_Pb"]
    mean_y = params["MEAN_Y_Pb"]

    sigma_x_vtx = params["SIGMA_X_VTX"]
    sigma_y_vtx = params["SIGMA_Y_VTX"]
    mean_x_vtx = params["MEAN_X_VTX"]
    mean_y_vtx = params["MEAN_Y_VTX"]
    MIN_CL_BEAM = params["MIN_CL_BEAM"]

    x0_min = -1000
    x0_max = 1000
    y0_min = -1000
    y0_max = 1000

    ntracks = [0,0,0,0,0,0]
    ntracks_sel = [0,0,0,0,0,0]


    eta_list = []
    eta_list_sel = []
    #Vertex QA plots
    res_tracks_pb = TH2D("res_tracks_pb","track position - Pb position at z = z_target; x_{trk}-x_{Pb} (mm); y_{trk}-y_{Pb} (mm)",2000,-5,5,2000,-5,5)
    res_tracks_pv = [
                        TH2D("res_tracks_pv","track position - PV position at z = z_target; x_{trk}-x_{V} (mm); y_{trk}-y_{V} (mm)",2000,-5,5,2000,-5,5),
                        TH2D("res_tracks_pv_ncl_2","track position - PV position at z = z_target; x_{trk}-x_{V} (mm); y_{trk}-y_{V} (mm)",2000,-5,5,2000,-5,5),
                        TH2D("res_tracks_pv_ncl_3","track position - PV position at z = z_target; x_{trk}-x_{V} (mm); y_{trk}-y_{V} (mm)",2000,-5,5,2000,-5,5),
                        TH2D("res_tracks_pv_ncl_4","track position - PV position at z = z_target; x_{trk}-x_{V} (mm); y_{trk}-y_{V} (mm)",2000,-5,5,2000,-5,5),
                        TH2D("res_tracks_pv_ncl_5","track position - PV position at z = z_target; x_{trk}-x_{V} (mm); y_{trk}-y_{V} (mm)",2000,-5,5,2000,-5,5),
                        TH2D("res_tracks_pv_ncl_6","track position - PV position at z = z_target; x_{trk}-x_{V} (mm); y_{trk}-y_{V} (mm)",2000,-5,5,2000,-5,5)
                    ]
    hEta = TH1D("hEta",";#eta;",100,0,10)
    
    hMultiplicity = [
                TH1D("hMultiplicity",";track multiplicity;",500,0,499.5),
                TH1D("hMultiplicity_ncl_2",";track multiplicity;",500,0,499.5),
                TH1D("hMultiplicity_ncl_3",";track multiplicity;",500,0,499.5),
                TH1D("hMultiplicity_ncl_4",";track multiplicity;",500,0,499.5),
                TH1D("hMultiplicity_ncl_5",";track multiplicity;",500,0,499.5),
                TH1D("hMultiplicity_ncl_6",";track multiplicity;",500,0,499.5)
            ]
    
    hMultiplicitySel = [
                TH1D("hMultiplicitySel",";track multiplicity;",500,0,499.5),
                TH1D("hMultiplicitySel_ncl_2",";track multiplicity;",500,0,499.5),
                TH1D("hMultiplicitySel_ncl_3",";track multiplicity;",500,0,499.5),
                TH1D("hMultiplicitySel_ncl_4",";track multiplicity;",500,0,499.5),
                TH1D("hMultiplicitySel_ncl_5",";track multiplicity;",500,0,499.5),
                TH1D("hMultiplicitySel_ncl_6",";track multiplicity;",500,0,499.5)
            ]
    hEtaSel = [
                TH1D("hEtaSel",";#eta;",100,0,10),
                TH1D("hEtaSel_ncl_2",";#eta;",100,0,10),
                TH1D("hEtaSel_ncl_3",";#eta;",100,0,10),
                TH1D("hEtaSel_ncl_4",";#eta;",100,0,10),
                TH1D("hEtaSel_ncl_5",";#eta;",100,0,10),
                TH1D("hEtaSel_ncl_6",";#eta;",100,0,10)
            ]

    hDCA =  [
                TH1D("hDCA",";DCA (mm);",1000,-10,10),
                TH1D("hDCA_ncl_2",";DCA (mm);",1000,-10,10),
                TH1D("hDCA_ncl_3",";DCA (mm);",1000,-1,1),
                TH1D("hDCA_ncl_4",";DCA (mm);",1000,-1,1),
                TH1D("hDCA_ncl_5",";DCA (mm);",1000,-1,1),
                TH1D("hDCA_ncl_6",";DCA (mm);",1000,-1,1)
            ]

    hDCAz = [
                TH1D("hDCAz",";DCA_{z} (mm);",1000,-1,1),
                TH1D("hDCAz_ncl_2",";DCA_{z} (mm);",1000,-1,1),
                TH1D("hDCAz_ncl_3",";DCA_{z} (mm);",1000,-0.1,0.1),
                TH1D("hDCAz_ncl_4",";DCA_{z} (mm);",1000,-0.1,0.1),
                TH1D("hDCAz_ncl_5",";DCA_{z} (mm);",1000,-0.1,0.1),
                TH1D("hDCAz_ncl_6",";DCA_{z} (mm);",1000,-0.1,0.1)
            ]
    hDCAxy = [
                TH1D("hDCAxy",";DCA_{xy} (mm);",1000,-10,10),
                TH1D("hDCAxy_ncl_2",";DCA_{xy} (mm);",1000,-10,10),
                TH1D("hDCAxy_ncl_3",";DCA_{xy} (mm);",1000,-1,1),
                TH1D("hDCAxy_ncl_4",";DCA_{xy} (mm);",1000,-1,1),
                TH1D("hDCAxy_ncl_5",";DCA_{xy} (mm);",1000,-1,1),
                TH1D("hDCAxy_ncl_6",";DCA_{xy} (mm);",1000,-1,1)
            ]
    
    hPbCluSizes = TH2D("hPbCluSizes",";cluster_size_1;cluster_size_2",200,-0.5,199.5,200,-0.5,199.5)
    hPbCluSizes_sel = TH2D("hPbCluSizes_sel",";cluster_size_1;cluster_size_2",200,-0.5,199.5,200,-0.5,199.5)
    hPbCluSize_1 = TH1D("hPbCluSize_ALPIDE0",";cluster_size;",200,-0.5,199.5)
    hPbCluSize_2 = TH1D("hPbCluSize_ALPIDE1",";cluster_size;",200,-0.5,199.5)
    hPbCluSize_1_sel = TH1D("hPbCluSize_ALPIDE0_sel",";cluster_size;",200,-0.5,199.5)
    hPbCluSize_2_sel = TH1D("hPbCluSize_ALPIDE1_sel",";cluster_size;",200,-0.5,199.5)
    hNcontr = TH1D("hNcontr",";number of contributors;",100,0,99.5)
    hNcontrVsZ = TH2D("hNcontrVsZ",";number of contributors;z_{v} (mm]);",100,0,99.5,1000,-2000,1000)
    hDeltaPbPVX = TH1D("hDeltaPbPVX",";x_{Pb}-x_{v} (mm);",1000,-0.5,0.5)
    hDeltaPbPVY = TH1D("hDeltaPbPVY",";y_{Pb}-y_{v} (mm);",1000,-0.5,0.5)
    hBeamVx = TH1D("hBeamVx","Beam direction x;beam v_{x};",1000,-0.1,0.1)
    hBeamVy = TH1D("hBeamVy","Beam direction y;beam v_{y};",1000,-0.1,0.1)
    hVertexX = TH1D("hVertexX",";x_{v} (mm);",1000,-10,10)
    hVertexY = TH1D("hVertexY",";y_{v} (mm);",1000,-10,10)
    hVertexZ = TH1D("hVertexZ",";z_{v} (mm);",100000,-2000,1000)
    hBeamVx_sel = TH1D("hBeamVx_sel","Beam direction x;beam v_{x};",1000,-0.1,0.1)
    hBeamVy_sel = TH1D("hBeamVy_sel","Beam direction y;beam v_{y};",1000,-0.1,0.1)
    hVertexX_sel = TH1D("hVertexX_sel",";x_{v} (mm);",1000,-10,10)
    hVertexY_sel = TH1D("hVertexY_sel",";y_{v} (mm);",1000,-10,10)
    hVertexZ_sel = TH1D("hVertexZ_sel",";z_{v} (mm);",100000,-2000,1000)
    hVertexZHighMult = TH1D("hVertexZHighMult",";z_{v} (mm);",100000,-2000,1000)
    file_ev = TFile(path, "r")
    hist_ev = file_ev.Get("EventLoaderEUDAQ2/ALPIDE_2/hPixelMultiplicityPerCorryEvent")
    nev_tot = hist_ev.GetEntries()
    #Vertexing
    df_vtx = pd.DataFrame(columns=df_tracks.columns.tolist())
    primary_vertexes = []
    # position (xv,yv,zv), nummber of contributors
    n_contributors = 0
    #loop over the index
    skip_event = False
    output = TFile("../results/output_na60plus_nsigma"+str(nsigma)+"_"+suffix+".root","recreate")

    # Record the start time
    start_time = time.time()
    if run_vertexing:
        for index, rows in df_tracks.iterrows():                   
            df_vtx = []

            pb_vx = df_beam.at[df_beam.index[index],"beam.vx"]
            pb_vy = df_beam.at[df_beam.index[index],"beam.vy"]
            pb_px = df_beam.at[df_beam.index[index],"beam.px"]
            pb_py = df_beam.at[df_beam.index[index],"beam.py"]
            
            hBeamVx.Fill(pb_vx)
            hBeamVy.Fill(pb_vy)

            n_contributors = 0
            x0 = pb_vx*z_target+pb_px
            y0 = pb_vy*z_target+pb_py
            if update:
                vertexer.define_PV_selection(res_tracks_pb, config, "Pb") 
            loading_bar(index+1, nev_tot)
            for vx,vy,px,py,cl1,cl2,cl3,cl4,cl5,cl6 in zip(rows["tracks.vx"],rows["tracks.vy"],rows["tracks.px"],rows["tracks.py"],rows["tracks.size1"],rows["tracks.size2"],rows["tracks.size3"],rows["tracks.size4"],rows["tracks.size5"],rows["tracks.size6"]):
                
                x = vx*z_target+px
                y = vy*z_target+py
                dx = ((x-x0)-mean_x)/sigma_x
                dy = ((y-y0)-mean_y)/sigma_y
                res_tracks_pb.Fill(x-x0,y-y0)

                cl_list = [cl1,cl2,cl3,cl4,cl5,cl6]
                n_clusters = count_clusters(cl_list)

                if dx**2+dy**2 < nsigma**2 and n_clusters > 2:#selecting the tracks 5 sigma from the primary vertex
                    n_contributors += 1
                    #row_to_add = df_tracks.iloc[0]  # Select the first row from df1
                    df_vtx.append([vx,vy,px,py,cl1,cl2,cl3,cl4,cl5,cl6])
        
            vertex = vertexer.fit_primary_vertex(df_vtx, df_beam.iloc[index], z_target)
            primary_vertexes.append([vertex[0],vertex[1],vertex[2],n_contributors])
            hDeltaPbPVX.Fill(x0-vertex[0])
            hDeltaPbPVY.Fill(y0-vertex[1])
            hVertexX.Fill(vertex[0])
            hVertexY.Fill(vertex[1])
            if n_contributors > 0:
                hVertexZ.Fill(vertex[2])
            if n_contributors > 20:
                hVertexZHighMult.Fill(vertex[2])
            hNcontr.Fill(n_contributors)
            hNcontrVsZ.Fill(n_contributors, vertex[2])

        np.save('primary_vertexes_'+suffix+'.npy', primary_vertexes)
        
        subdir_vtx = output.mkdir("vertex-qa")
        subdir_vtx.cd()
        res_tracks_pb.Write()
        hDeltaPbPVX.Write()
        hDeltaPbPVY.Write()
        hBeamVx.Write()
        hBeamVy.Write()
        hVertexX.Write()
        hVertexY.Write()
        hVertexZ.Write()
        hVertexZHighMult.Write()
        hNcontr.Write()
        hNcontrVsZ.Write()
    else:
        primary_vertexes = np.load('primary_vertexes_'+suffix+'.npy')

    if run_analysis:
        #analysis
        for index, rows in df_tracks.iterrows(): 

            skip_event = False

            max_clsiz = max(df_beam.at[df_beam.index[index],"beam.size1"], df_beam.at[df_beam.index[index],"beam.size2"])
            
            hPbCluSize_1.Fill(df_beam.at[df_beam.index[index],"beam.size1"])
            hPbCluSize_2.Fill(df_beam.at[df_beam.index[index],"beam.size2"])
            hPbCluSizes.Fill(df_beam.at[df_beam.index[index],"beam.size1"],df_beam.at[df_beam.index[index],"beam.size2"])
            if max_clsiz < MIN_CL_BEAM:
                skip_event = True
            else:
                hPbCluSizes_sel.Fill(df_beam.at[df_beam.index[index],"beam.size1"],df_beam.at[df_beam.index[index],"beam.size2"])
                hPbCluSize_1_sel.Fill(df_beam.at[df_beam.index[index],"beam.size1"])
                hPbCluSize_2_sel.Fill(df_beam.at[df_beam.index[index],"beam.size2"])

            x0, y0, z0, n_con = primary_vertexes[index]

            if z0 < 80 or z0 > 120 or n_con == 0:
                skip_event = True
            if x0 < x0_min or x0 > x0_max:
                skip_event = True
            if y0 < y0_min or y0 > y0_max:
                skip_event = True
            if not skip_event:
                hBeamVx_sel.Fill(df_beam.at[df_beam.index[index],"beam.vx"])
                hBeamVy_sel.Fill(df_beam.at[df_beam.index[index],"beam.vy"])
                hVertexX_sel.Fill(x0)
                hVertexY_sel.Fill(y0)
                hVertexZ_sel.Fill(z0)

            for cl in range(6):
                hMultiplicity[cl].Fill(ntracks[cl])
                hMultiplicitySel[cl].Fill(ntracks_sel[cl])

            ntracks = [0,0,0,0,0,0]
            ntracks_sel = [0,0,0,0,0,0]
                
            loading_bar(index+1, nev_tot, process="Analysis")
            if skip_event:
                continue
            for vx,vy,px,py,cl1,cl2,cl3,cl4,cl5,cl6 in zip(rows["tracks.vx"],rows["tracks.vy"],rows["tracks.px"],rows["tracks.py"],rows["tracks.size1"],rows["tracks.size2"],rows["tracks.size3"],rows["tracks.size4"],rows["tracks.size5"],rows["tracks.size6"]):
                
                x = vx*z_target+px
                y = vy*z_target+py
                res_tracks_pb.Fill(x-x0,y-y0)

                cl_list = [cl1,cl2,cl3,cl4,cl5,cl6]
                n_clusters = count_clusters(cl_list)
                
                ntracks[0] += 1
                ntracks[n_clusters-1] += 1

                dx = ((x-x0)-mean_x_vtx)/sigma_x_vtx
                dy = ((y-y0)-mean_y_vtx)/sigma_y_vtx

                eta = math.atanh(1./math.sqrt(vx**2+vy**2+1))

                if dx**2+dy**2 < 3**2: #selecting the tracks 5 sigma from the primary vertex
                    ntracks_sel[0] += 1
                    ntracks_sel[n_clusters-1] += 1
                    hEtaSel[0].Fill(eta)
                    hEtaSel[n_clusters-1].Fill(eta)
                hEta.Fill(eta)
                distances = vertexer.compute_dca([x0,y0,z0],[vx,vy,px,py])
                dca = math.sqrt(distances[0]**2+distances[1]**2+distances[2]**2)
                dcaz = distances[2]
                dcaxy = math.sqrt(distances[0]**2+distances[1]**2)

                res_tracks_pv[0].Fill(x-x0, y-y0)
                res_tracks_pv[n_clusters-1].Fill(x-x0, y-y0)

                hDCA[0].Fill(dca)
                hDCA[n_clusters-1].Fill(dca)

                hDCAz[0].Fill(dcaz)
                hDCAz[n_clusters-1].Fill(dcaz)

                hDCAxy[0].Fill(dcaxy)
                hDCAxy[n_clusters-1].Fill(dcaxy)

        if update:        
            vertexer.define_PV_selection(res_tracks_pv, config, "VTX")

        subdir_mult = output.mkdir("multiplicity")

        subdir_pb = subdir_mult.mkdir("beam-qa")
        subdir_pb.cd()
        hPbCluSize_1.Write()
        hPbCluSize_2.Write()
        hPbCluSizes.Write()
        hPbCluSizes_sel.Write()
        hPbCluSize_1_sel.Write()
        hPbCluSize_2_sel.Write()
        hBeamVx_sel.Write()
        hBeamVy_sel.Write()
        hVertexX_sel.Write()
        hVertexY_sel.Write()
        hVertexZ_sel.Write()

        subdir_dca = subdir_mult.mkdir("dca")
        subdir_dca.cd()
        for i in range(6):
            hDCA[i].Write()
        for i in range(6):
            hDCAz[i].Write()
        for i in range(6):
            hDCAxy[i].Write()
        subdir_mult.cd()
        for i in range(6):
            hEtaSel[i].Write()  
        for i in range(6):      
            res_tracks_pv[i].Write()
        for i in range(6):
            hMultiplicity[i].Write()
        for i in range(6):
            hMultiplicitySel[i].Write()
        hEta.Write()
        
        output.Close()

    # Record the end time
    end_time = time.time()
    # Calculate the elapsed time
    elapsed_time = end_time - start_time
    print(f"\nExcution time: {elapsed_time} seconds")

def plot_results(path_list, data_list, config_list):
    legend_multiplicity = TLegend(0.7,0.7,0.9,0.9)
    color_list = [
                   ROOT.kRed,
                   ROOT.kBlue,
                   ROOT.kGreen,
                   ROOT.kOrange,
                   ROOT.kBlack,
                      
                ]
    cv_multiplicity = TCanvas("cv_multiplicity", "cv_multiplicity")

    # Create a list to store the cloned histograms
    cloned_histograms = []

    for path, data, config, color in zip(path_list, data_list, config_list, color_list):
        fIn = TFile(path)
        hMulti = fIn.Get("multiplicity/hMultiplicitySel")
        file_ev = TFile(data, "r")
        hist_ev = file_ev.Get("EventLoaderEUDAQ2/ALPIDE_2/hPixelMultiplicityPerCorryEvent")
        nev_tot = hist_ev.GetEntries()

        with open(os.path.expandvars(config), 'r') as stream:
            try:
                params = yaml.full_load(stream)
            except yaml.YAMLError as exc:
                print(exc)

        TARGET = params["TARGET"]
        MASS_NUMBER = params["MASS_NUMBER"]

        # Create a copy of the histogram
        hMultiClone = hMulti.Clone("hMultiClone_" + TARGET)
        hMultiClone.Scale(1.0 / (MASS_NUMBER * nev_tot))
        hMultiClone.SetLineColor(color)
        hMultiClone.GetYaxis().SetTitle("event/(#trigger x A)")
        hMultiClone.GetXaxis().SetRangeUser(0, 150)

        cv_multiplicity.cd()
        hMultiClone.Draw("same")
        cloned_histograms.append(hMultiClone)  # Store the clone in the list
        legend_multiplicity.AddEntry(hMultiClone, TARGET, "lep")

    legend_multiplicity.Draw()
    cv_multiplicity.SetLogy()
    cv_multiplicity.SaveAs("multiplicity.png")





def main():
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument("-r","--read_tree", help="Read data and produce output histograms", action="store_true")
    parser.add_argument("-u","--update", help="Update sigmas for PV selection", action="store_true")
    parser.add_argument("-v","--vertexing", help="Run the vertexing", action="store_true")
    parser.add_argument("-a","--analysis", help="Run the analysis", action="store_true")
    parser.add_argument("-p","--plot_results", help="Produce plot from the analysis", action="store_true")
    parser.add_argument("config", help="Path to the YAML configuration file")
    args = parser.parse_args()


    if args.read_tree:
        suffix = "422231250231017232859_Be_planes_2_Unique_False_dEta0.06_pruned"
        path = "/home/giacomo/NA60+_test_beam_data/NA60Plus_test_beam/analysis_"+suffix+".root"

        read_tree(path= path,compute_sigmas=True,config=args.config,nsigma=3, suffix=suffix, update=args.update, run_vertexing=args.vertexing, run_analysis=args.analysis)

    if args.plot_results:
        path_Ag = "/home/giacomo/its-corryvreckan-tools/output/2023-10_SPS/Target/analysis_422213212231017213217_Ag_planes_6.root"
        path_Air = "/home/giacomo/its-corryvreckan-tools/output/2023-10_SPS/Target/analysis_422191053231017193343_Pb_planes_6.root"
        path_Be = "/home/giacomo/its-corryvreckan-tools/output/2023-10_SPS/Target/analysis_422231249231017231254_Be_planes_6.root"
        path_Pb = "/home/giacomo/NA60+_test_beam_data/NA60Plus_test_beam/analysis_423133844231018133850_False_planes_6.root"
        path_S = "/home/giacomo/NA60+_test_beam_data/NA60Plus_test_beam/analysis_423140916231018140921_S_planes_6.root"
        path_list = [
                        "/home/giacomo/NA60+_test_beam_data/NA60Plus_test_beam/results/output_na60plus_nsigma3_Sulfur.root",
                        "/home/giacomo/NA60+_test_beam_data/NA60Plus_test_beam/results/output_na60plus_nsigma3_ARg.root"
                    ]
            
        config_list = [
                        "/home/giacomo/NA60+_test_beam_data/NA60Plus_test_beam/macros/configs/linear_setup_S.yaml",
                        "/home/giacomo/NA60+_test_beam_data/NA60Plus_test_beam/macros/configs/linear_setup_Ag.yaml"
                    ]
        
        data_list = [
                        "/home/giacomo/NA60+_test_beam_data/NA60Plus_test_beam/analysis_423140916231018140921_S_planes_6.root",
                        "/home/giacomo/NA60+_test_beam_data/NA60Plus_test_beam/analysis_422213213231017215647_Ag_planes_6.root"
                    ]
            
        plot_results(path_list, data_list, config_list)

main()
