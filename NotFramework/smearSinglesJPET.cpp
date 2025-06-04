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

void smearSinglesJPET(TString fileName)
{
    TPRegexp regexp("\\bSingles(\\w+)");
    TString ss = fileName(regexp);
    ss.Remove(0,7);
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
    
    double time;
    int rsectorID, crystalID;
    float posX, posY, posZ;
    
    input_tree->SetBranchAddress("time",&time);
    input_tree->SetBranchAddress("rsectorID",&rsectorID);
    input_tree->SetBranchAddress("crystalID",&crystalID);
    input_tree->SetBranchAddress("posX",&posX);
    input_tree->SetBranchAddress("posY",&posY);
    input_tree->SetBranchAddress("posZ",&posZ);
    
    input_tree->BuildIndex("0", "1e15 * time");
    
    TTreeIndex *index = (TTreeIndex*)input_tree->GetTreeIndex();
    
    TString outName;
    outName.Form("sorted_Singles%i.root",fileNumber);    
    TFile f2(outName,"RECREATE");
    
    TTree *tsorted = (TTree*)input_tree->CloneTree(0);
    
    fstream centers("modular_centers.txt", ios::in);
    if(!centers.is_open())
    {
        cout << "File not opened..." << endl;
    }
    std::vector <double> vx;
    std::vector <double> vy;
    while(true)
    {
        double tx, ty;
        centers >> tx >> ty;
        if(centers.eof())break;
        vx.push_back(tx);
        vy.push_back(ty);
    }
    centers.close();
    
    gRandom->SetSeed(0);
    
    for( long int i = 0; i < index->GetN(); i++ )
    {
        long int local = input_tree->LoadTree( index->GetIndex()[i] );
        input_tree->GetEntry(local);
        int id = crystalID + 13*rsectorID;
        posX = (float)vx[id];
        posY = (float)vy[id];
        // --- IN CASE OF Z SMEARING ---
        /*
        float z = gRandom->Gaus( posZ , 15. );
        posZ = z;
        */
        // --- ---
        time = time;
        if(posZ<250. && posZ>-250.)tsorted->Fill();
    }
    tsorted->Write();
}
