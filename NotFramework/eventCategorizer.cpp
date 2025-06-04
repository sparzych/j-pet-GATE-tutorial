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

#include "TVector3.h"

using namespace std;

bool noiseCut(float);
double getDistance(TLorentzVector v1, TLorentzVector v2);
bool scatterTest(TLorentzVector v1, TLorentzVector v2);
std::pair<double, double> scatterDiff(TLorentzVector v1, TLorentzVector v2);
TVector3 calculateAnnihilationPoint(TLorentzVector v1, TLorentzVector v2);
double calculateAngle(TLorentzVector v1, TLorentzVector v2);

void eventCategorizer(TString fileName)
{
    TPRegexp regexp("\\bTimeWindows(\\w+)");
    TString ss = fileName(regexp);
    ss.Remove(0,11);
    int fileNumber = ss.Atoi();
    
    TFile *input = new TFile(fileName);
    TTree *input_tree = (TTree*)input->Get("TimeWindows");
    
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
    input_tree->SetBranchStatus("multiplicity",1);
    
    std::vector<int> *event = 0;
    std::vector<int> *photon = 0;
    std::vector<int> *source = 0;
    std::vector<int> *rsector = 0;
    std::vector<int> *crystal = 0;
    std::vector<int> *cC = 0;
    std::vector<int> *cP = 0;
    std::vector<int> *rC = 0;
    std::vector<int> *rP = 0;
    std::vector<float> *energy = 0;
    std::vector<float> *x = 0;
    std::vector<float> *y = 0;
    std::vector<float> *z = 0;
    std::vector<float> *sx = 0;
    std::vector<float> *sy = 0;
    std::vector<float> *sz = 0;
    std::vector<double> *time = 0;
    int multiplicity;
    
    input_tree->SetBranchAddress("time",&time);
    input_tree->SetBranchAddress("rsectorID",&rsector);
    input_tree->SetBranchAddress("crystalID",&crystal);
    input_tree->SetBranchAddress("eventID",&event);
    input_tree->SetBranchAddress("photonID",&photon);
    input_tree->SetBranchAddress("source",&source);
    input_tree->SetBranchAddress("comptonPhantom",&cP);
    input_tree->SetBranchAddress("comptonCrystal",&cC);
    input_tree->SetBranchAddress("RayleighPhantom",&rP);
    input_tree->SetBranchAddress("RayleighCrystal",&rC);
    input_tree->SetBranchAddress("energy",&energy);
    input_tree->SetBranchAddress("posX",&x);
    input_tree->SetBranchAddress("posY",&y);
    input_tree->SetBranchAddress("posZ",&z);
    input_tree->SetBranchAddress("sX",&sx);
    input_tree->SetBranchAddress("sY",&sy);
    input_tree->SetBranchAddress("sZ",&sz);
    input_tree->SetBranchAddress("multiplicity",&multiplicity); 
    
    TH1D *scat_test_time = new TH1D("scat_test_time","scat_test_time",1000,-5000,5000);
    
    for(long int i = 0; i < input_tree->GetEntries(); i++)
    {
        input_tree->GetEntry(i);
        std::vector<int> goodHits = {};
        
        for(int j = 0; j < multiplicity; j++)
        {
            if( noiseCut( energy->at(j) ) )
            {
                goodHits.push_back(j);
            }
        }
        
        if( goodHits.size() == 2 )
        {
            int id1 = goodHits.at(0);
            int id2 = goodHits.at(1);
                                    
            TLorentzVector v1( x->at(id1), y->at(id1), z->at(id1), time->at(id1) );
            TLorentzVector v2( x->at(id2), y->at(id2), z->at(id2), time->at(id2) );
            
            float t1x = x->at(id1) - 0.;
            float t1y = y->at(id1) - 0.;
            float t1z = z->at(id1) - 0.;
            TVector3 track1( t1x, t1y, t1z );
            
            float t2x = x->at(id2) - x->at(id1);
            float t2y = y->at(id2) - y->at(id1);
            float t2z = z->at(id2) - z->at(id1);
            TVector3 track2( t2x, t2y, t2z );
            
            double testTime, testDist;
            std::tie(testTime, testDist) = scatterDiff(v1, v2);
            scat_test_time->Fill( testTime );
        }
    }
       
    TString outName;
    outName.Form("histos_%i.root",fileNumber);
    TFile histos(outName,"RECREATE");
    scat_test_time->Write();
    histos.Close();
}

bool noiseCut(float e)
{
    if(e > 0.03)
        return true;
    return false;
}

double getDistance(TLorentzVector v1, TLorentzVector v2)
{
    float d = sqrt( pow(v1.X() - v2.X(), 2.) + pow(v1.Y() - v2.Y(), 2.) + pow(v1.Z() - v2.Z(), 2.) );
    return (double)d;
}

double calculateAngle(TLorentzVector v1, TLorentzVector v2)
{
    TVector3 d1( v1.X(), v1.Y(), v1.Z() );
    TVector3 d2( v2.X(), v2.Y(), v2.Z() );
    return TMath::RadToDeg() * ( d1.Angle(d2) );
}

bool scatterTest(TLorentzVector v1, TLorentzVector v2)
{
    double c = 299792458000.;   // mm/s
    double dist = getDistance(v1, v2);  // [mm]
    double timeDiff = v2.T() - v1.T();  // [s]
    double testDist = dist - timeDiff*c;    // [mm]
    double testTime = timeDiff - dist/c;    // [s]
    
    testTime *= pow(10.,12);    // [ps]
    
    if( testTime < -1000. )
        return true;
    return false;
}

std::pair<double, double> scatterDiff(TLorentzVector v1, TLorentzVector v2)
{
    double c = 299792458000.;   // mm/s
    double dist = getDistance(v1, v2);  // [mm]
    double timeDiff = v2.T() - v1.T();  // [s]
    double testDist = dist - timeDiff*c;    // [mm]
    double testTime = timeDiff - dist/c;    // [s]
    
    testTime *= pow(10.,12);    // [ps]
    testDist /= 10.;     // [cm]
    
    return std::make_pair(testTime, testDist);
}

TVector3 calculateAnnihilationPoint(TLorentzVector v1, TLorentzVector v2)
{
    TVector3 hit1(v1.X(), v1.Y(), v1.Z());
    TVector3 hit2(v2.X(), v2.Y(), v2.Z());
    
    TVector3 middleOfLOR = 0.5 * (hit1 + hit2);
    TVector3 versorOnLOR = (hit2 - hit1).Unit();
    
    double c = 299792458000.;   // mm/s
    double tof = v1.T() - v2.T();
    double shift = 0.5 * tof * c;
    
    TVector3 annihilationPoint;
    annihilationPoint.SetX( middleOfLOR.X() + shift * versorOnLOR.X() );
    annihilationPoint.SetY( middleOfLOR.Y() + shift * versorOnLOR.Y() );
    annihilationPoint.SetZ( middleOfLOR.Z() + shift * versorOnLOR.Z() );
    
    return annihilationPoint;
}
