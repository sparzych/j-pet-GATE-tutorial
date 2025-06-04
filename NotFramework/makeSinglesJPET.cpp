#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"

using namespace std;

void makeSinglesJPET(TString fileName)
{
    TPRegexp regexp("\\bresults(\\w+)");
    TString ss = fileName(regexp);
    ss.Remove(0,7);
    int fileNumber = ss.Atoi();
    
    TChain *s = new TChain("Singles");
    s->Add(fileName);
    
    s->SetBranchStatus("*",0);
    s->SetBranchStatus("eventID",1);
    s->SetBranchStatus("rsectorID",1);
    s->SetBranchStatus("crystalID",1);
    s->SetBranchStatus("layerID",1);
    s->SetBranchStatus("moduleID",1);
    s->SetBranchStatus("submoduleID",1);
    s->SetBranchStatus("comptonPhantom",1);
    s->SetBranchStatus("RayleighPhantom",1);
    s->SetBranchStatus("comptonCrystal",1);
    s->SetBranchStatus("RayleighCrystal",1);
    s->SetBranchStatus("energy",1);
    s->SetBranchStatus("time",1);
    s->SetBranchStatus("globalPosX",1);
    s->SetBranchStatus("globalPosY",1);
    s->SetBranchStatus("globalPosZ",1);
    s->SetBranchStatus("sourcePosX",1);
    s->SetBranchStatus("sourcePosY",1);
    s->SetBranchStatus("sourcePosZ",1);
    
    int eventID, rsectorID, crystalID, layerID, comptonPhantom, RayleighPhantom, moduleID, submoduleID, comptonCrystal, RayleighCrystal;
    float energy, posX, posY, posZ, sX, sY, sZ;
    double time;
    
    s->SetBranchAddress("eventID",&eventID);
    s->SetBranchAddress("rsectorID",&rsectorID);
    s->SetBranchAddress("crystalID",&crystalID);
    s->SetBranchAddress("layerID",&layerID);
    s->SetBranchAddress("moduleID",&moduleID);
    s->SetBranchAddress("submoduleID",&submoduleID);
    s->SetBranchAddress("comptonPhantom",&comptonPhantom);
    s->SetBranchAddress("RayleighPhantom",&RayleighPhantom);
    s->SetBranchAddress("comptonCrystal",&comptonCrystal);
    s->SetBranchAddress("RayleighCrystal",&RayleighCrystal);
    s->SetBranchAddress("energy",&energy);
    s->SetBranchAddress("time",&time);
    s->SetBranchAddress("globalPosX",&posX);
    s->SetBranchAddress("globalPosY",&posY);
    s->SetBranchAddress("globalPosZ",&posZ);
    s->SetBranchAddress("sourcePosX",&sX);
    s->SetBranchAddress("sourcePosY",&sY);
    s->SetBranchAddress("sourcePosZ",&sZ);
    
    TChain *h = new TChain("Hits");
    h->Add(fileName);
    
    h->SetBranchStatus("*",0);
    h->SetBranchStatus("eventID",1);
    h->SetBranchStatus("rsectorID",1);
    h->SetBranchStatus("crystalID",1);
    h->SetBranchStatus("layerID",1);
    h->SetBranchStatus("moduleID",1);
    h->SetBranchStatus("submoduleID",1);
    h->SetBranchStatus("photonID",1);
    h->SetBranchStatus("gammaType",1);
    h->SetBranchStatus("sourceType",1);
    
    int HeventID, HrsectorID, HcrystalID, HlayerID, HphotonID, HmoduleID, HsubmoduleID, HgammaType, HsourceType;
    
    h->SetBranchAddress("eventID",&HeventID);
    h->SetBranchAddress("rsectorID",&HrsectorID);
    h->SetBranchAddress("crystalID",&HcrystalID);
    h->SetBranchAddress("layerID",&HlayerID);
    h->SetBranchAddress("moduleID",&HmoduleID);
    h->SetBranchAddress("submoduleID",&HsubmoduleID);
    h->SetBranchAddress("photonID",&HphotonID);
    h->SetBranchAddress("gammaType",&HgammaType);
    h->SetBranchAddress("sourceType",&HsourceType);
    
    TString treeName;
    treeName.Form("Singles%i.root",fileNumber);
    TFile *file = new TFile(treeName,"RECREATE");
    int SphotonID, ScomptonPhantom, SRayleighPhantom, SeventID, SrsectorID, ScomptonCrystal, SRayleighCrystal, Ssource, ScrystalID;
    float Senergy, SposX, SposY, SposZ, SsX, SsY, SsZ;
    double Stime;
    TTree *tree = new TTree("Singles","Singles");
    tree->Branch("eventID",&SeventID,"eventID/I");
    tree->Branch("photonID",&SphotonID,"photonID/I");
    tree->Branch("source",&Ssource,"source/I");
    tree->Branch("rsectorID",&SrsectorID,"rsectorID/I");
    tree->Branch("crystalID",&ScrystalID,"crystalID/I");
    tree->Branch("time",&Stime,"time/D");
    tree->Branch("energy",&Senergy,"energy/F");
    tree->Branch("posX",&SposX,"posX/F");
    tree->Branch("posY",&SposY,"posY/F");
    tree->Branch("posZ",&SposZ,"posZ/F");
    tree->Branch("sX",&SsX,"sX/F");
    tree->Branch("sY",&SsY,"sY/F");
    tree->Branch("sZ",&SsZ,"sZ/F");
    tree->Branch("comptonPhantom",&ScomptonPhantom,"comptonPhantom/I");
    tree->Branch("RayleighPhantom",&SRayleighPhantom,"RayleighPhantom/I");
    tree->Branch("comptonCrystal",&ScomptonCrystal,"comptonCrystal/I");
    tree->Branch("RayleighCrystal",&SRayleighCrystal,"RayleighCrystal/I");
    
    long int start = 0;
    
    long int diff = 0;
    
    for(int i = 0; i < s->GetEntries(); i++)
    {
        s->GetEntry(i);
        int flag = 0;
        int flag_photon = 0;
        int previous = -1;
        int save = 0;
        int ph_ID = -1;
        int source_ID = -1;
        
        for(int j = start; j< h->GetEntries(); j++)
        {
            h->GetEntry(j);
            if(HeventID == eventID)
            {
                if(flag==0)
                {
                    start = j;
                }
                flag = 1;
                
                if(HrsectorID == rsectorID && HcrystalID == crystalID && HlayerID == layerID && HmoduleID == moduleID/* && HsubmoduleID == submoduleID*/ )//&& HgammaType != 0)
                {
                    if(flag_photon == 0)
                    {
                        previous = HphotonID;
                        flag_photon = 1;
                    }
                    
                    if( (HphotonID != previous) && (HphotonID - previous != 3) )
                    {
                        cout << "Not matching photons in " << HeventID << endl;
                        diff++;
                        cout << previous << "\t" << HphotonID << endl;
                        save = 2;
                    }
                    else
                    {
                        if(save != 2)
                        {
                            source_ID = HsourceType;
                            ph_ID = HphotonID;
                            save = 1;
                        }
                    }
                }
            }
            else
            {
                if(flag == 1)
                {
                    flag = 0;
                    flag_photon = 0;
                    break;
                }
            }
        }
        if(save == 1)
        {
            Ssource = source_ID;
            SphotonID = ph_ID;
            Senergy = energy;
            ScomptonPhantom = comptonPhantom;
            SRayleighPhantom = RayleighPhantom;
            ScomptonCrystal = comptonCrystal;
            SRayleighCrystal = RayleighCrystal;
            
            // --- IN CASE OF TIME SMEARING WITH DEPENDECE ON TOT ---
            /*
            double A0 = 5.648;
            double A1 = 7.433;
            double A2 = 0.994;
            double tot = A0 - A1 * pow(A2, (energy*1000.));
            tot *= 1000.;
            double t_res = 5.41642*pow(10.,5) /sqrt(tot + 1.46080*pow(10.,4) ) - 3.62140*pow(10.,3);
            Stime = gRandom->Gaus( time , t_res/pow(10.,12) );
            */
            // --- ---
            Stime = time;
            
            SeventID = eventID;
            SrsectorID = rsectorID;
            ScrystalID = crystalID;
            SposX = posX;
            SposY = posY;
            SposZ = posZ;
            SsX = sX;
            SsY = sY;
            SsZ = sZ;
            if(energy>0.001)tree->Fill();
            save = 0;
            ph_ID = -1;
        }
    }
    tree->Write();
    file->Close();
    cout << "Sum of strange: " << diff << endl;
}
