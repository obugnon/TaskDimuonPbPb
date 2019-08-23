/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* AliAnaysisTaskMyTask
 *
 * empty task which can serve as a starting point for building an analysis
 * as an example, one histogram is filled
 */

// include root libraries
#include <iostream>
#include "TChain.h"
#include "TH1F.h"
#include "TList.h"
#include "TChain.h"
#include "TMath.h"
#include "THnSparse.h"
// include AliRoot Libraries
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliVEvent.h"
#include "AliAODInputHandler.h"
#include "AliMuonTrackCuts.h"
#include "AliVTrack.h"
#include "AliAnalysisMuonUtility.h"
#include "AliMultSelection.h"
#include "TaskDimuonPbPb.h"

class TaskDimuonPbPb;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(TaskDimuonPbPb) // classimp: necessary for root

TaskDimuonPbPb::TaskDimuonPbPb() : AliAnalysisTaskSE(), 
    fAODEvent(0),
    fVEvent(0),
    fMuonTrackCuts(0),
    fTriggerClass(0),
    fFirstRun(0),
    fLastRun(0),
    fListEventHistos(0), 
    fListSingleMuonHistos(0),
    fListDiMuonHistos(0),
    fHistoTotalEventsPerRun(0),
    fHistoPSEventsPerRun(0),
    fHistoEventsBeforePSPerRun(0),
    fHistoCMULEventsInCINT7(0),
    fHistoCMSLEventsInCINT7(0),
    fHistoCMULEventsInCMSL(0),
    fHistoNumberMuonsCuts(0),
    fHistoDiMuonOS(0),
    fHistoDiMuonLS(0),
    fHistoSingleMuon(0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
TaskDimuonPbPb::TaskDimuonPbPb(const char* name,int firstRun, int lastRun, UInt_t triggerClass) : AliAnalysisTaskSE(name),
    fAODEvent(0),
    fVEvent(0),
    fMuonTrackCuts(0),
    fTriggerClass(triggerClass),
    fFirstRun(firstRun),
    fLastRun(lastRun),
    fListEventHistos(0), 
    fListSingleMuonHistos(0),
    fListDiMuonHistos(0),
    fHistoTotalEventsPerRun(0),
    fHistoPSEventsPerRun(0),
    fHistoEventsBeforePSPerRun(0),
    fHistoCMULEventsInCINT7(0),
    fHistoCMSLEventsInCINT7(0),
    fHistoCMULEventsInCMSL(0),
    fHistoNumberMuonsCuts(0),
    fHistoDiMuonOS(0),
    fHistoDiMuonLS(0),
    fHistoSingleMuon(0)
{
    // constructor
    DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
                                        // this chain is created by the analysis manager, so no need to worry about it, 
                                        // it does its work automatically
    DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms 
    DefineOutput(2, TList::Class());    // you can add more output objects by calling DefineOutput(2, classname::Class())
    DefineOutput(3, TList::Class());    // if you add more output objects, make sure to call PostData for all of them, and to
                                        // make changes to your AddTask macro!
}
//_____________________________________________________________________________
TaskDimuonPbPb::~TaskDimuonPbPb()
{
    // destructor
    if(fListEventHistos && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
        delete fListEventHistos;     // at the end of your task, it is deleted from memory by calling this function
    }
    if(fListSingleMuonHistos && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
        delete fListSingleMuonHistos;
    }
    if(fListDiMuonHistos && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
        delete fListDiMuonHistos;
    }
}
//_____________________________________________________________________________
void TaskDimuonPbPb::NotifyRun()
{
  /// Set run number for cuts
  if ( fMuonTrackCuts ) fMuonTrackCuts->SetRun(fInputHandler);
}
//_____________________________________________________________________________
void TaskDimuonPbPb::UserCreateOutputObjects()
{
    // create output objects
    //
    // this function is called ONCE at the start of your analysis (RUNTIME)
    // here you ceate the histograms that you want to use 
    //
    // the histograms are in this case added to a tlist, this list is in the end saved
    // to an output file

     //Event histograms
    fListEventHistos = new TList();
    fListEventHistos->SetOwner(kTRUE);

    fHistoTotalEventsPerRun = new TH1I("fHistoTotalEventsPerRun","",fLastRun - fFirstRun,fFirstRun,fLastRun);
    fHistoTotalEventsPerRun->Sumw2();
    fHistoTotalEventsPerRun->GetXaxis()->SetTitle("Run Number");
    fHistoTotalEventsPerRun->GetYaxis()->SetTitle("# Total Events");
    fListEventHistos->Add(fHistoTotalEventsPerRun);
    
    fHistoPSEventsPerRun = new TH1I("fHistoPSEventsPerRun","",fLastRun - fFirstRun,fFirstRun,fLastRun);
    fHistoPSEventsPerRun->Sumw2();
    fHistoPSEventsPerRun->GetXaxis()->SetTitle("Run Number");
    fHistoPSEventsPerRun->GetYaxis()->SetTitle("# PS Events");
    fListEventHistos->Add(fHistoPSEventsPerRun);


    fHistoEventsBeforePSPerRun = new TH1I("fHistoEventsBeforePSPerRun","",fLastRun - fFirstRun, fFirstRun, fLastRun);
    fHistoEventsBeforePSPerRun->Sumw2();
    fHistoEventsBeforePSPerRun->GetXaxis()->SetTitle("Run Number");
    fHistoEventsBeforePSPerRun->GetYaxis()->SetTitle("# Events before physic selection");
    fListEventHistos->Add(fHistoEventsBeforePSPerRun);

    if(fTriggerClass == AliVEvent::kINT7inMUON)
    {
        fHistoCMULEventsInCINT7 = new TH1I("fHistoCMULEventsInCINT7","",fLastRun - fFirstRun,fFirstRun,fLastRun);
        fHistoCMULEventsInCINT7->Sumw2();
        fHistoCMULEventsInCINT7->GetXaxis()->SetTitle("Run Number");
        fHistoCMULEventsInCINT7->GetYaxis()->SetTitle("# 0MUL inputs in CINT7 events");
        fListEventHistos->Add(fHistoCMULEventsInCINT7);
    
        fHistoCMSLEventsInCINT7 = new TH1I("fHistoCMSLEventsInCINT7","",fLastRun - fFirstRun,fFirstRun,fLastRun);
        fHistoCMSLEventsInCINT7->Sumw2();
        fHistoCMSLEventsInCINT7->GetXaxis()->SetTitle("Run Number");
        fHistoCMSLEventsInCINT7->GetYaxis()->SetTitle("# 0MSL inputs in CINT7 events");
        fListEventHistos->Add(fHistoCMSLEventsInCINT7);
    }


    if(fTriggerClass == AliVEvent::kMuonSingleLowPt7)
    {
        fHistoCMULEventsInCMSL = new TH1I("fHistoCMULEventsInCMSL","",fLastRun - fFirstRun,fFirstRun,fLastRun);
        fHistoCMULEventsInCMSL->Sumw2();
        fHistoCMULEventsInCMSL->GetXaxis()->SetTitle("Run Number");
        fHistoCMULEventsInCMSL->GetYaxis()->SetTitle("# 0MUL inputs in CMSL events");
        fListEventHistos->Add(fHistoCMULEventsInCMSL);
    }

    fListSingleMuonHistos = new TList();
    fListSingleMuonHistos->SetOwner(kTRUE);
    fListDiMuonHistos = new TList();
    fListDiMuonHistos->SetOwner(kTRUE);

    if(fTriggerClass == AliVEvent::kMuonUnlikeLowPt7)
    {
        //SingleMuon histograms
        Int_t nbinsSingleMuon[5]={1000,60,100,100, 20}; //pT, Eta, Theta, Phi, cent
        Double_t xminSingleMuon[5]={0,-5,0.75*TMath::Pi(),-TMath::Pi(), 0}, xmaxSingleMuon[5]={100,-2,1.25*TMath::Pi(),TMath::Pi(),100};
        fHistoSingleMuon = new THnSparseD("fHistoSingleMuon","",5, nbinsSingleMuon,xminSingleMuon,xmaxSingleMuon, 1024*16);
        fHistoSingleMuon->Sumw2();
        fHistoSingleMuon->GetAxis(0)->SetTitle("p_{T} GeV/c");
        fHistoSingleMuon->GetAxis(1)->SetTitle("#eta");
        fHistoSingleMuon->GetAxis(2)->SetTitle("#theta");
        fHistoSingleMuon->GetAxis(3)->SetTitle("#phi");
        fHistoSingleMuon->GetAxis(4)->SetTitle("Centrality of the event %");
        fListSingleMuonHistos->Add(fHistoSingleMuon);

        fHistoNumberMuonsCuts = new TH1F("fHistoNumberMuonsCuts","",2,0,2);
        fHistoNumberMuonsCuts->Sumw2();
        fHistoNumberMuonsCuts->GetXaxis()->SetTitle("bin 1: without cut, bin 2: with cuts");
        fHistoNumberMuonsCuts->GetYaxis()->SetTitle("# single muons");
        fListSingleMuonHistos->Add(fHistoNumberMuonsCuts);


        //DiMuon histograms
        Int_t nbinsDiMuon[4]={1000,400,60,20}; //Mmumu, pT, y, centrality
        Double_t xminDiMuon[4]={0,0,-5,0}, xmaxDiMuon[4]={20,20,-2,100};
        fHistoDiMuonOS = new THnSparseD("fHistoDiMuonOS","",4,nbinsDiMuon,xminDiMuon,xmaxDiMuon, 1024*16);
        fHistoDiMuonOS->Sumw2();
        fHistoDiMuonOS->GetAxis(0)->SetTitle("M_{#mu#mu} GeV/c^{2}");
        fHistoDiMuonOS->GetAxis(1)->SetTitle("p_{T} GeV/c");
        fHistoDiMuonOS->GetAxis(2)->SetTitle("y");
        fHistoDiMuonOS->GetAxis(3)->SetTitle("Centrality of the event %");
        fListDiMuonHistos->Add(fHistoDiMuonOS);

        fHistoDiMuonLS = new THnSparseD("fHistoDiMuonLS","",4,nbinsDiMuon,xminDiMuon,xmaxDiMuon, 1024*16);
        fHistoDiMuonLS->Sumw2();
        fHistoDiMuonLS->GetAxis(0)->SetTitle("M_{#mu#mu} GeV/c^{2}");
        fHistoDiMuonLS->GetAxis(1)->SetTitle("p_{T} GeV/c");
        fHistoDiMuonLS->GetAxis(2)->SetTitle("y");
        fHistoDiMuonLS->GetAxis(3)->SetTitle("Centrality of the event %");
        fListDiMuonHistos->Add(fHistoDiMuonLS);


        //The muon muonTrackCuts can be defined here. Hiwever it is better to defien it outside (in addTaskDimuonPPB.C). To be fixed
        fMuonTrackCuts = new AliMuonTrackCuts("StandardMuonTrackCuts","StandardMuonTrackCuts");
        fMuonTrackCuts->SetAllowDefaultParams(kTRUE);
        fMuonTrackCuts->SetFilterMask (AliMuonTrackCuts::kMuEta | AliMuonTrackCuts::kMuThetaAbs  | AliMuonTrackCuts::kMuMatchLpt | AliMuonTrackCuts::kMuPdca);//Set the cuts to be used for the muon selections. See all the available cuts in AliMuonTrackCuts.h
    }

  //This is needed to save the outputs.
  PostData(1, fListEventHistos);
  PostData(2, fListSingleMuonHistos);
  PostData(3, fListDiMuonHistos);

}
//_____________________________________________________________________________
void TaskDimuonPbPb::UserExec(Option_t *)
{
    // user exec this function is called once for each event
    // the manager will take care of reading the events from file, and with the static function InputEvent() you 
    // have access to the current event. 
    // once you return from the UserExec function, the manager will retrieve the next event from the chain
    fAODEvent = dynamic_cast<AliAODEvent*>(InputEvent());    
    if(!fAODEvent) {
        AliError("ERROR: Could not retrieve AOD event !!");
        return;
    }
    int runNumber;
    fVEvent = static_cast<AliVEvent *>(InputEvent());
    runNumber = fAODEvent->GetRunNumber();
    
    AliMultSelection *multSelection = (AliMultSelection * ) fAODEvent->FindListObject("MultSelection");
    Double_t centralityFromV0 = multSelection->GetMultiplicityPercentile("V0M", false);
    if (centralityFromV0 > 90.) return;

    // Event Histos   
    fHistoTotalEventsPerRun->Fill(runNumber);
    TString strFiredTriggers = fAODEvent->GetFiredTriggerClasses();

    if(fTriggerClass == AliVEvent::kINT7inMUON)
    {
        if(strFiredTriggers.Contains("CINT7-B-NOPF-MUFAST")) fHistoEventsBeforePSPerRun->Fill(runNumber);
    }

    if(fTriggerClass == AliVEvent::kMuonUnlikeLowPt7)
    {
        if(strFiredTriggers.Contains("CMUL7-B-NOPF-MUFAST")) fHistoEventsBeforePSPerRun->Fill(runNumber);
    }

    if(fTriggerClass == AliVEvent::kMuonSingleLowPt7)
    {
        if(strFiredTriggers.Contains("CMSL7-B-NOPF-MUFAST")) fHistoEventsBeforePSPerRun->Fill(runNumber);
    }
    
                                
    //Physics Selection
    UInt_t IsSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected());
    if(IsSelected & fTriggerClass)
    {
        fHistoPSEventsPerRun->Fill(runNumber);
        UInt_t L0InputMUL = 1<<20;
        UInt_t L0InputV0A = 1<<0;
        UInt_t L0InputV0C = 1<<1;
        UInt_t L0InputMSL = 1<<17;

        AliAODHeader *aodHeader = NULL;
        aodHeader = (AliAODHeader*)fAODEvent->GetHeader();
        UInt_t L0TriggerInputs = aodHeader->GetL0TriggerInputs();

        if(fTriggerClass == AliVEvent::kINT7inMUON)
        {
            if( (L0TriggerInputs & L0InputV0A) && (L0TriggerInputs & L0InputV0C) )//MB condition (should be the same as the condition for isSelectedCINTEvent, just to make sure)
            {
                if(L0TriggerInputs & L0InputMUL){fHistoCMULEventsInCINT7->Fill(runNumber);}
                if(L0TriggerInputs & L0InputMSL){fHistoCMSLEventsInCINT7->Fill(runNumber);}
            }
        }

        if(fTriggerClass == AliVEvent::kMuonSingleLowPt7)
        {
            if( (L0TriggerInputs & L0InputV0A) && (L0TriggerInputs & L0InputV0C) && (L0TriggerInputs & L0InputMSL))//CMSL condition
            {
                if(L0TriggerInputs & L0InputMUL){fHistoCMULEventsInCMSL->Fill(runNumber);}
            }
        }

        //Fill Single Muon and Dimuon properties histograms
        if(fTriggerClass == AliVEvent::kMuonUnlikeLowPt7)
        {   
            TLorentzVector lvMuon1, lvMuon2, lvDiMuon;
            AliVParticle* muonTrack1;
            AliVParticle* muonTrack2;

            Float_t muonMass2 = AliAnalysisMuonUtility::MuonMass2(); // the PDG rest mass of the muon (constante, used for getting the kinematics) en GeV
            int numberOfTracks = AliAnalysisMuonUtility::GetNTracks(fVEvent); // get the number of muon tracks in the event
            for(Int_t iMuon1 = 0 ; iMuon1 < numberOfTracks ; iMuon1++) // loop ove rall these tracks
            {
                muonTrack1 = AliAnalysisMuonUtility::GetTrack(iMuon1,fVEvent);
                if( !muonTrack1 ) { AliError(Form("ERROR: Could not retrieve AOD or ESD track %d", iMuon1)); continue;}

                //Number muons before/after muon cuts
                fHistoNumberMuonsCuts->Fill(0.5);
                if ( ! fMuonTrackCuts->IsSelected(muonTrack1) ) continue;//include cuts on pDCA, Eta, Rabs
                fHistoNumberMuonsCuts->Fill(1.5);

                //single muon properties 
                Float_t energy1 = TMath::Sqrt(muonTrack1->P()*muonTrack1->P() + muonMass2);
                lvMuon1.SetPxPyPzE(muonTrack1->Px(),muonTrack1->Py(),muonTrack1->Pz(),energy1); //def 4-vect muon1
                Short_t muonCharge1 = muonTrack1->Charge();
                
                Double_t propertiesSingleMuon[5]={lvMuon1.Pt(),lvMuon1.Eta(),lvMuon1.Theta(),lvMuon1.Phi(),centralityFromV0};
                fHistoSingleMuon->Fill(propertiesSingleMuon,1);
                
                for (Int_t iMuon2 = iMuon1+1; iMuon2 < numberOfTracks; iMuon2++)
                {
                    muonTrack2 = AliAnalysisMuonUtility::GetTrack(iMuon2,fVEvent);
                    if ( !muonTrack2 ) {AliError(Form("ERROR: Could not retrieve AOD or ESD track %d", iMuon2)); continue;}
                    if ( ! fMuonTrackCuts->IsSelected(muonTrack2) ) continue;//include cuts on pDCA, Eta, Rabs

                    Float_t energy2 = TMath::Sqrt(muonTrack2->P()*muonTrack2->P() + muonMass2);
                    lvMuon2.SetPxPyPzE(muonTrack2->Px(),muonTrack2->Py(),muonTrack2->Pz(),energy2); //def 4-vect muon1
                    Short_t muonCharge2 = muonTrack2->Charge();

                    //dimuon properties
                    lvDiMuon = lvMuon1+lvMuon2;

                    Double_t propertiesDiMuon[4]={};
                    propertiesDiMuon[0]=lvDiMuon.M();
                    propertiesDiMuon[1]=lvDiMuon.Pt();
                    propertiesDiMuon[2]=lvDiMuon.Rapidity();
                    propertiesDiMuon[3]=centralityFromV0;

                    if (muonCharge1 != muonCharge2){ fHistoDiMuonOS->Fill(propertiesDiMuon,1); }
                    else if (muonCharge1 == muonCharge2){ fHistoDiMuonLS->Fill(propertiesDiMuon,1); }

                }
                      
            }
        }
    }

    PostData(1, fListEventHistos);
    PostData(2, fListSingleMuonHistos);
    PostData(3, fListDiMuonHistos);
}
//_____________________________________________________________________________
void TaskDimuonPbPb::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
