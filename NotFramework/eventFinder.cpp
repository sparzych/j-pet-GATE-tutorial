#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TH1F.h"
#include "TBenchmark.h"
#include "TRandom.h"
#include "TSystem.h"
#include "TString.h"

using namespace std;

void eventFinder(TString fileName)
{
    TPRegexp regexp("\\bsorted_Singles(\\w+)");
    TString ss = fileName(regexp);
    ss.Remove(0,14);
    int fileNumber = ss.Atoi();
    
    TFile *input = new TFile(fileName);
    TTree *input_tree = (TTree*)input->Get("Singles");
    
    input_tree->SetBranchStatus("*",0);
    input_tree->SetBranchStatus("time",1);
    input_tree->SetBranchStatus("rsectorID",1);
    input_tree->SetBranchStatus("crystalID",1);
    input_tree->SetBranchStatus("eventID",1);
    input_tree->SetBranchStatus("photonID",1);
    input_tree->SetBranchStatus("source",1);
    input_tree->SetBranchStatus("comptonPhantom",1);
    input_tree->SetBranchStatus("comptonCrystal",1);
    input_tree->SetBranchStatus("RayleighPhantom",1);
    input_tree->SetBranchStatus("RayleighCrystal",1);
    input_tree->SetBranchStatus("energy",1);
    input_tree->SetBranchStatus("posX",1);
    input_tree->SetBranchStatus("posY",1);
    input_tree->SetBranchStatus("posZ",1);
    input_tree->SetBranchStatus("sX",1);
    input_tree->SetBranchStatus("sY",1);
    input_tree->SetBranchStatus("sZ",1);
    
    double Itime;
    int rsectorID, eventID, comptonPhantom, comptonCrystal, RayleighPhantom, RayleighCrystal, photonID, Isource, crystalID;
    float Ienergy, posX, posY, posZ, sX, sY, sZ;
    
    input_tree->SetBranchAddress("time",&Itime);
    input_tree->SetBranchAddress("rsectorID",&rsectorID);
    input_tree->SetBranchAddress("crystalID",&crystalID);
    input_tree->SetBranchAddress("eventID",&eventID);
    input_tree->SetBranchAddress("photonID",&photonID);
    input_tree->SetBranchAddress("source",&Isource);
    input_tree->SetBranchAddress("comptonPhantom",&comptonPhantom);
    input_tree->SetBranchAddress("comptonCrystal",&comptonCrystal);
    input_tree->SetBranchAddress("RayleighPhantom",&RayleighPhantom);
    input_tree->SetBranchAddress("RayleighCrystal",&RayleighCrystal);
    input_tree->SetBranchAddress("energy",&Ienergy);
    input_tree->SetBranchAddress("posX",&posX);
    input_tree->SetBranchAddress("posY",&posY);
    input_tree->SetBranchAddress("posZ",&posZ);
    input_tree->SetBranchAddress("sX",&sX);
    input_tree->SetBranchAddress("sY",&sY);
    input_tree->SetBranchAddress("sZ",&sZ);
    
    TString outName;
    outName.Form("TimeWindows%i.root",fileNumber);
    TFile *output = new TFile(outName,"RECREATE");
    TTree *output_tree = new TTree("TimeWindows","Tree with time windows");
    
    std::vector<int> event;
    std::vector<int> photon;
    std::vector<int> source;
    std::vector<int> rsector;
    std::vector<int> crystal;
    std::vector<int> cC;
    std::vector<int> cP;
    std::vector<int> rC;
    std::vector<int> rP;
    std::vector<float> energy;
    std::vector<float> x;
    std::vector<float> y;
    std::vector<float> z;
    std::vector<float> sx;
    std::vector<float> sy;
    std::vector<float> sz;
    std::vector<double> time;
    int m;
    
    output_tree->Branch("eventID",&event);
    output_tree->Branch("photonID",&photon);
    output_tree->Branch("source",&source);
    output_tree->Branch("rsectorID",&rsector);
    output_tree->Branch("crystalID",&crystal);
    output_tree->Branch("comptonPhantom",&cP);
    output_tree->Branch("RayleighPhantom",&rP);
    output_tree->Branch("comptonCrystal",&cC);
    output_tree->Branch("RayleighCrystal",&rC);
    output_tree->Branch("energy",&energy);
    output_tree->Branch("posX",&x);
    output_tree->Branch("posY",&y);
    output_tree->Branch("posZ",&z);
    output_tree->Branch("sX",&sx);
    output_tree->Branch("sY",&sy);
    output_tree->Branch("sZ",&sz);
    output_tree->Branch("time",&time);
    output_tree->Branch("multiplicity",&m);
    
    double window = 3.0 * pow(10.,-9.);
    
    input_tree->GetEntry(0);
    double start = Itime;
        
    for(long int i = 0; i < input_tree->GetEntries(); i++)
    {
        input_tree->GetEntry(i);
//         cout << std::setprecision(8) << Itime << endl;
        if(Itime > start + window)
        {
            int window_steps = (int)((Itime - start)/window);
            start += window_steps*window;
            
            m = event.size();
            if(m!=0)
            {
                output_tree->Fill();
            }
                        
            event.clear();
            photon.clear();
            source.clear();
            rsector.clear();
            crystal.clear();
            cP.clear();
            rP.clear();
            cC.clear();
            rC.clear();
            energy.clear();
            x.clear();
            y.clear();
            z.clear();
            sx.clear();
            sy.clear();
            sz.clear();
            time.clear();
        }
    
        event.emplace_back(eventID);
        photon.emplace_back(photonID);
        source.emplace_back(Isource);
        rsector.emplace_back(rsectorID);
        crystal.emplace_back(crystalID);
        cP.emplace_back(comptonPhantom);
        rP.emplace_back(RayleighPhantom);
        cC.emplace_back(comptonCrystal);
        rC.emplace_back(RayleighCrystal);
        energy.emplace_back(Ienergy);
        x.emplace_back(posX);
        y.emplace_back(posY);
        z.emplace_back(posZ);
        sx.emplace_back(sX);
        sy.emplace_back(sY);
        sz.emplace_back(sZ);
        time.emplace_back(Itime);
    }
    
    output->Write();
    delete input;
    delete output;
}
